#include <chrono>
#include <filesystem>
#include <sstream>
#include <map>
#include <set>

#include "src/builder/build.cpp"
#include "src/dictionary.cpp"
#include "src/info.cpp"
#include "include/streaming_query.hpp"

using namespace sshash;

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename Kmer, typename Offsets>
void check(std::filesystem::path const& directory, std::vector<std::string> const& strings,
           uint64_t k, uint64_t m, bool check_stream) {
    dictionary<Kmer, Offsets> dict;
    dict.call_minimizer_positions(0, [](uint64_t) {
        throw std::runtime_error("empty dictionary emitted a minimizer position");
    });
    auto input = directory / "input.fa";
    {
        std::ofstream out(input);
        for (auto const& string : strings) out << ">test\n" << string << '\n';
    }
    build_configuration config;
    config.k = k;
    config.m = m;
    config.verbose = false;
    config.allow_reverse_complement_duplicates = true;
    config.tmp_dirname = directory.string();
    dict.build(input.string(), config);

    for (unsigned int byte = 128; byte < 256; ++byte) {
        require(!Kmer::is_valid(static_cast<char>(byte)), "non-ASCII byte accepted as alphabet character");
        auto query = strings.front().substr(0, k);
        query.front() = static_cast<char>(byte);
        streaming_query<decltype(dict)> stream(&dict);
        require(stream.lookup(query.data()).kmer_id == constants::invalid_uint64,
                "streaming lookup accepted non-ASCII input");
    }

    auto verify = [&] {
        std::map<uint64_t, std::set<uint64_t>> expected_positions;
        uint64_t base = 0;
        uint64_t string_id = 0;
        for (auto const& string : strings) {
            auto [begin, end] = dict.string_offsets(string_id++);
            require(begin == base && end == base + string.size(), "incorrect string boundaries");
            for (uint64_t i = 0; i + k <= string.size(); ++i) {
                auto word = util::string_to_uint_kmer<Kmer>(string.data() + i, k);
                auto mini = util::compute_minimizer(word, k, m, dict.hasher());
                expected_positions[mini.minimizer].insert(base + i + mini.pos_in_kmer);
            }
            base += string.size();
        }
        // Include absent canonical keys, not just the minimizers seen above.
        for (uint64_t mini = 0; mini < (uint64_t(1) << (Kmer::bits_per_char * m)); ++mini) {
            if (util::canonical_mmer<Kmer>(mini, m) != mini) continue;
            std::set<uint64_t> observed;
            dict.call_minimizer_positions(mini, [&](uint64_t offset) {
                require(observed.insert(offset).second, "minimizer enumeration repeated a position");
            });
            auto expected = expected_positions.find(mini);
            require(expected == expected_positions.end() ? observed.empty() : observed == expected->second,
                    "minimizer enumeration lost or added indexed positions");
        }
        uint64_t expected_id = 0;
        for (auto const& string : strings) {
            streaming_query<decltype(dict)> stream(&dict);
            for (uint64_t i = 0; i + k <= string.size(); ++i, ++expected_id) {
                auto result = dict.lookup(string.data() + i, false);
                if (result.kmer_id != expected_id) {
                    std::cerr << "k=" << k << " m=" << m << " word=" << string.substr(i, k)
                              << " expected=" << expected_id << " got=" << result.kmer_id << '\n';
                }
                require(result.kmer_id == expected_id, "forward lookup lost a stored orientation");
                require(result.kmer_orientation == constants::forward_orientation,
                        "stored kmer must have forward orientation");
                if (check_stream) {
                    require(equal_lookup_result(result, stream.lookup(string.data() + i)),
                            "streaming lookup differs from direct lookup");
                }
            }
        }
    };
    verify();
    std::stringstream bytes;
    essentials::generic_saver saver(bytes);
    saver.visit(dict);
    essentials::generic_loader loader(bytes);
    loader.visit(dict);
    verify();
}

template <typename Offsets>
void check_dna(std::filesystem::path const& directory) {
    using Kmer = dna_uint_kmer_t<uint64_t>;
    check<Kmer, Offsets>(directory, {"ACGT"}, 3, 1, false);
    // Endpoints at powers of two need an extra bit; exercise neighboring sizes
    // and multiple records, including after serialization.
    for (auto sequence : {"AAA", "AAAC", "AAACC", "AAACCGT", "AAACCGTA", "AAACCGTAT"}) {
        check<Kmer, Offsets>(directory, {sequence}, 3, 1, false);
    }
    check<Kmer, Offsets>(directory, {"AAAC", "CCGT"}, 3, 1, false);
    // One record per kmer forces heavy buckets; both DNA orientations share
    // canonical skew-index keys but must retain their distinct forward IDs.
    std::vector<std::string> kmers;
    for (uint64_t value = 0; value < (1ULL << 10); ++value) {
        std::string word(5, 'A');
        for (uint64_t i = 0; i < 5; ++i) word[i] = "ACGT"[(value >> (2 * i)) & 3];
        kmers.push_back(word);
    }
    check<Kmer, Offsets>(directory, kmers, 5, 1, false);

    dictionary<Kmer, Offsets> dict;
    auto input = directory / "negative.fa";
    { std::ofstream out(input); out << ">test\nAAAC\n"; }
    build_configuration config;
    config.k = 3;
    config.m = 1;
    config.verbose = false;
    config.tmp_dirname = directory.string();
    dict.build(input.string(), config);
    require(dict.lookup("TTT", false).kmer_id == constants::invalid_uint64,
            "forward-only lookup accepted reverse complement");
    require(dict.lookup("TTT").kmer_id == 0, "reverse complement lookup failed");
}

int main() {
    auto directory = std::filesystem::temp_directory_path() /
        ("sshash-lookup-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directory(directory);
    try {
        check_dna<decoded_offsets>(directory);
        check_dna<encoded_offsets>(directory);
        const std::vector<std::string> protein = {"MKWVTFISRLDPQEAHNYCGMKWVTFISRLDPQEAHNYCG"};
        // Distinct windows; the long kmer also exercises a multiword iterator.
        check<aa_uint_kmer_t<uint64_t>, decoded_offsets>(directory,
            {"MKWVTFISRLDPQEAHNYCG"}, 3, 1, true);
        check<aa_uint_kmer_t<__uint128_t>, encoded_offsets>(directory, protein, 21, 3, true);
        // Dense protein buckets require more control bits than sparse buckets.
        std::vector<std::string> protein_kmers;
        const std::string alphabet = "ACDEFGHIKLMNPQRSTVWY";
        for (char a : alphabet) for (char b : alphabet) for (char c : alphabet) {
            protein_kmers.push_back(std::string{a, b, c});
        }
        check<aa_uint_kmer_t<uint64_t>, decoded_offsets>(directory, protein_kmers, 3, 1, true);
        std::filesystem::remove_all(directory);
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        std::filesystem::remove_all(directory);
        return 1;
    }
}
