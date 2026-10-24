#include "include/kmer.hpp"
#include <array>
#include <cstring>
#include "include/util.hpp"
#include "external/gz/zip_stream.hpp"

#if defined(__AVX2__)
#include <immintrin.h>
#include <x86intrin.h>
#endif

namespace sshash {

namespace util {

#if defined(__AVX2__)
/*
    Pack 32 two-bit values, one per byte, into a single 64-bit word.
*/
inline uint64_t pack2bits(__m256i values) {
    // collect bit-0 plane
    __m256i bit0 = _mm256_slli_epi16(values, 7);
    uint32_t mask0 = _mm256_movemask_epi8(bit0);

    // collect bit-1 plane
    __m256i bit1 = _mm256_slli_epi16(values, 6);
    uint32_t mask1 = _mm256_movemask_epi8(bit1);

    // interleave into the 64-bit result
    uint64_t even = _pdep_u64(mask0, 0x5555555555555555ULL);  // 010101...
    uint64_t odd = _pdep_u64(mask1, 0xAAAAAAAAAAAAAAAAULL);   // 101010...
    return even | odd;
}
#endif

template <typename Kmer>
inline uint64_t pack_word(char const* sequence) {
#if defined(__AVX2__)
    if constexpr (Kmer::bits_per_char == 2) {
        uint8_t codes[32];
        for (uint64_t j = 0; j != 32; ++j) {
            codes[j] = Kmer::char_to_uint(sequence[j]);
        }
        return pack2bits(_mm256_loadu_si256(reinterpret_cast<__m256i const*>(codes)));
    }
#endif
    uint64_t word = 0;
    for (uint64_t j = 0; j != 64 / Kmer::bits_per_char; ++j) {
        word |= uint64_t(Kmer::char_to_uint(sequence[j])) << (j * Kmer::bits_per_char);
    }
    return word;
}

template <typename Kmer>
struct batch_mapper {
    static uint64_t convert(char c) { return Kmer::char_to_uint(c); }
};
// Keep scalar lookup table-based; this form lets the bulk mapping loop vectorize.
template <typename Storage>
struct batch_mapper<aa_uint_kmer_t<Storage>> {
    static uint64_t convert(char c) {
        const unsigned char code = (static_cast<unsigned char>(c) | 0x20) - 'a';
        return code < 26 ? code + 2 * (code == 23) - (code >= 24) : uint64_t(-1);
    }
};
template <size_t N>
inline void write_words(uint64_t* destination, unsigned shift, std::array<uint64_t, N> const& words) {
    if (shift == 0) {
        std::copy(words.begin(), words.end(), destination);
    } else {
        destination[0] |= words[0] << shift;
        for (size_t j = 1; j < N; ++j)
            destination[j] = (words[j] << shift) | (words[j - 1] >> (64 - shift));
        destination[N] = words[N - 1] >> (64 - shift);
    }
}
template <typename Kmer, size_t N>
inline auto pack_batch(char const* sequence) {
    constexpr unsigned B = Kmer::bits_per_char;
    static_assert(N % 64 == 0);
    uint8_t codes[N];
    for (size_t j = 0; j != N; ++j) codes[j] = batch_mapper<Kmer>::convert(sequence[j]);
    uint64_t groups[N / 8];
    using lanes [[gnu::vector_size(32)]] = uint64_t;
    lanes blocks[N / 32];
    std::memcpy(blocks, codes, N);
    for (auto& words : blocks) {
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
        for (unsigned i = 0; i < 4; ++i) words[i] = __builtin_bswap64(words[i]);
#endif
        auto merge = [&](unsigned group, uint64_t mask) {
            words = (words & mask) | ((words >> (group * (8 - B))) & (mask << (group * B)));
        };
        // Merge adjacent groups of one, two, and four byte lanes.
        merge(1, ((uint64_t(1) << B) - 1) * 0x0001000100010001ULL);
        merge(2, ((uint64_t(1) << (2 * B)) - 1) * 0x0000000100000001ULL);
        merge(4, (uint64_t(1) << (4 * B)) - 1);
    }
    std::memcpy(groups, blocks, N);
    std::array<uint64_t, N * B / 64> result{};
    // A 128-character batch produces at most 16 words; unroll to keep shifts constant.
#pragma GCC unroll 16
    for (size_t i = 0; i < result.size(); ++i) {
        const size_t first = 64 * i / (8 * B);
        result[i] = groups[first] >> (64 * i % (8 * B));
        for (size_t j = first + 1; j * (8 * B) < 64 * (i + 1); ++j)
            result[i] |= groups[j] << (j * (8 * B) - 64 * i);
    }
    return result;
}
// Keep the bulk loop out of encode_strings so short records retain its small loop.
template <typename Kmer>
[[gnu::noinline]] void append_record_batches(bits::bit_vector::builder& output, char const* sequence, uint64_t n) {
    constexpr uint64_t batch_chars = 128;
    constexpr uint64_t words_per_batch = batch_chars * Kmer::bits_per_char / 64;
    const uint64_t count = n / batch_chars;
    const uint64_t old_bits = output.num_bits();
    const unsigned shift = old_bits % 64;
    // Grow once per record, then write whole words, carrying across its bit offset.
    output.resize(old_bits + count * words_per_batch * 64);
    uint64_t* destination = output.data().data() + old_bits / 64;
    for (uint64_t j = 0; j != count; ++j) {
        auto words = pack_batch<Kmer, batch_chars>(sequence + j * batch_chars);
        write_words(destination + j * words_per_batch, shift, words);
    }
    uint64_t i = count * batch_chars;
    constexpr uint64_t chars_per_word = 64 / Kmer::bits_per_char;
    for (; i + chars_per_word <= n; i += chars_per_word)
        output.append_bits(pack_word<Kmer>(sequence + i), chars_per_word * Kmer::bits_per_char);
    for (; i < n; ++i) output.append_bits(Kmer::char_to_uint(sequence[i]), Kmer::bits_per_char);
}

}  // namespace util

template <typename Kmer, typename Offsets>
void dictionary_builder<Kmer, Offsets>::encode_strings(std::istream& is,
                                                       const input_file_t fmt)  //
{
    const uint64_t k = build_config.k;
    const uint64_t m = build_config.m;
    assert(k > 0 and k >= m);
    const uint64_t max_num_kmers_in_super_kmer = k - m + 1;

    if (max_num_kmers_in_super_kmer >= (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8))) {
        throw std::runtime_error(
            "max_num_kmers_in_super_kmer " + std::to_string(max_num_kmers_in_super_kmer) +
            " does not fit into " + std::to_string(sizeof(num_kmers_in_super_kmer_uint_type) * 8) +
            " bits");
    }

    {
        const uint64_t num_bits_for_strings = 8 * 8 * essentials::GB;  // reserve 8 GB of memory
        const uint64_t num_sequences = 100'000'000;
        strings_builder.reserve(num_bits_for_strings);
        strings_offsets_builder.reserve(num_sequences);
    }

    std::string sequence;
    uint64_t num_bases = 0;
    uint64_t max_len = 0;
    uint64_t seq_len = 0;
    weights_builder.init();

    /* intervals of weights */
    uint64_t weight_value = constants::invalid_uint64;
    uint64_t weight_length = 0;

    while (true)  //
    {
        if (fmt == input_file_t::cf_seg) {
            std::getline(is, sequence, '\t');  // skip until '\t' and consume it
        } else {
            assert(fmt == input_file_t::fasta);
            if (build_config.weighted) {     // parse header
                std::getline(is, sequence);  // header sequence
                if (sequence.empty()) break;

                /*
                    Heder format:
                    >[id] LN:i:[seq_len] ab:Z:[weight_seq]
                    where [weight_seq] is a space-separated sequence of integer counters
                    (the weights), whose length is equal to [seq_len]-k+1.
                    Example: '>12 LN:i:41 ab:Z:2 2 2 2 2 2 2 2 2 2 2'
                */

                expect(sequence[0], '>');
                uint64_t i = 0;
                i = sequence.find_first_of(' ', i);
                if (i == std::string::npos) throw parse_runtime_error();

                i += 1;
                expect(sequence[i + 0], 'L');
                expect(sequence[i + 1], 'N');
                expect(sequence[i + 2], ':');
                expect(sequence[i + 3], 'i');
                expect(sequence[i + 4], ':');
                i += 5;
                uint64_t j = sequence.find_first_of(' ', i);
                if (j == std::string::npos) throw parse_runtime_error();

                seq_len = std::strtoull(sequence.data() + i, nullptr, 10);
                i = j + 1;
                expect(sequence[i + 0], 'a');
                expect(sequence[i + 1], 'b');
                expect(sequence[i + 2], ':');
                expect(sequence[i + 3], 'Z');
                expect(sequence[i + 4], ':');
                i += 5;

                for (uint64_t j = 0; j != seq_len - k + 1; ++j) {
                    uint64_t weight = std::strtoull(sequence.data() + i, nullptr, 10);
                    i = sequence.find_first_of(' ', i) + 1;
                    weights_builder.eat(weight);
                    if (weight == weight_value) {
                        weight_length += 1;
                    } else {
                        if (weight_value != constants::invalid_uint64) {
                            weights_builder.push_weight_interval(weight_value, weight_length);
                        }
                        weight_value = weight;
                        weight_length = 1;
                    }
                }
            } else {
                // skip header sequence
                is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
        }

        std::getline(is, sequence);  // DNA sequence

        if (is.eof()) break;

        const uint64_t n = sequence.length();
        assert(n >= k);
        max_len = n > max_len ? n : max_len;
        assert(strings_builder.num_bits() % Kmer::bits_per_char == 0);
        strings_offsets_builder.push_back(strings_builder.num_bits() / Kmer::bits_per_char);
        num_kmers += n - k + 1;
        num_bases += n;

        if (build_config.weighted and seq_len != n) {
            std::cout << "ERROR: expected a sequence of length " << seq_len
                      << " but got one of length " << n << std::endl;
            throw std::runtime_error("file is malformed");
        }

        if (build_config.verbose and strings_offsets_builder.size() % 1'000'000 == 0) {
            std::cout << "read " << strings_offsets_builder.size() << " sequences, " << num_bases
                      << " bases, " << num_kmers << " kmers" << std::endl;
        }

        if constexpr (Kmer::bits_per_char != 2 && Kmer::bits_per_char <= 8) {
            if (n >= 128) {
                util::append_record_batches<Kmer>(strings_builder, sequence.data(), n);
                continue;
            }
        }
        uint64_t i = 0;
        constexpr uint64_t chars_per_word = 64 / Kmer::bits_per_char;
        constexpr uint64_t bits_per_word = chars_per_word * Kmer::bits_per_char;
        for (; i + chars_per_word <= n; i += chars_per_word) {
            strings_builder.append_bits(util::pack_word<Kmer>(sequence.data() + i), bits_per_word);
        }
        for (; i < n; ++i) {
            strings_builder.append_bits(Kmer::char_to_uint(sequence[i]), Kmer::bits_per_char);
        }
    }

    strings_offsets_builder.push_back(strings_builder.num_bits() / Kmer::bits_per_char);
    assert(strings_offsets_builder.front() == 0);
    assert(strings_offsets_builder.size() >= 2);

    /* Push a final sentinel (dummy) value to avoid bounds' checking in
       kmer_iterator::fill_buff(). */
    static_assert(Kmer::uint_kmer_bits % 64 == 0);
    for (int dummy_bits = Kmer::uint_kmer_bits; dummy_bits; dummy_bits -= 64) {
        strings_builder.append_bits(0, 64);
    }

    const uint64_t num_sequences = strings_offsets_builder.size() - 1;

    if (build_config.verbose) {
        std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases, "
                  << num_kmers << " kmers" << std::endl;
        std::cout << "num_kmers " << num_kmers << std::endl;
        std::cout << "cost: 2.0 + "
                  << static_cast<double>(Kmer::bits_per_char * num_sequences * (k - 1)) / num_kmers
                  << " [bits/kmer]" << std::endl;
    }

    /*
        The parameter m (minimizer length) should be at least
            ceil(log_s(N))+1
        where N is the number of nucleotides in the input and s is the alphabet size.
        We warn the user if the used m is less than this lower bound.
    */
    const uint64_t s = uint64_t(1) << Kmer::bits_per_char;
    const uint64_t lower_bound_on_m = std::ceil(std::log(num_bases) / std::log(s)) + 1;
    if (build_config.verbose and m < lower_bound_on_m) {
        std::cout << "\n--> WARNING: using minimizer length " << m
                  << " but it should be at least ceil(log_" << s << "(" << num_bases
                  << "))+1 = " << lower_bound_on_m << '\n'
                  << std::endl;
    }

    if (build_config.weighted) {
        weights_builder.push_weight_interval(weight_value, weight_length);
        weights_builder.finalize(num_kmers);
    }

    num_bits nb;
    nb.per_absolute_offset = std::ceil(std::log2(strings_offsets_builder.back()));
    nb.per_relative_offset = std::ceil(std::log2(max_len - m + 1));
    nb.per_string_id = std::ceil(std::log2(num_sequences));

    if (build_config.verbose) {
        std::cout << "max string length = " << max_len << std::endl;
        std::cout << "num bits per_absolute_offset = " << nb.per_absolute_offset << std::endl;
        std::cout << "num bits per_relative_offset = " << nb.per_relative_offset << std::endl;
        std::cout << "num bits per_string_id = " << nb.per_string_id << std::endl;
    }

    if (nb.per_string_id + nb.per_relative_offset > 64) {
        throw std::runtime_error("minimizer offset does not fit within 64 bits");
    }

    strings_offsets_builder.set_num_bits(nb);
}

template <typename Kmer, typename Offsets>
void dictionary_builder<Kmer, Offsets>::encode_strings(std::string const& filename)  //
{
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    if (build_config.verbose) std::cout << "reading file '" << filename << "'..." << std::endl;
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        if (util::ends_with(filename, ".cf_seg.gz")) {
            encode_strings(zis, input_file_t::cf_seg);
        } else {
            encode_strings(zis, input_file_t::fasta);
        }
    } else {
        if (util::ends_with(filename, ".cf_seg")) {
            encode_strings(is, input_file_t::cf_seg);
        } else {
            encode_strings(is, input_file_t::fasta);
        }
    }
    is.close();
}

}  // namespace sshash
