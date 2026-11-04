#include <chrono>
#include <filesystem>
#include <random>
#include <unordered_map>

#include "src/builder/build.cpp"
#include "src/dictionary.cpp"
#include "src/info.cpp"

using namespace sshash;

struct test_directory {
    std::filesystem::path path = std::filesystem::temp_directory_path() /
        ("sshash-heavy-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    test_directory() { std::filesystem::create_directory(path); }
    ~test_directory() {
        std::error_code error;
        std::filesystem::remove_all(path, error);
    }
};

template <typename Offsets>
void check(std::filesystem::path const& directory, std::vector<std::string> const& words) {
    auto input = directory / "input.fa";
    {
        std::ofstream out(input);
        out.exceptions(std::ios::failbit | std::ios::badbit);
        for (auto const& word : words) out << ">test\n" << word << '\n';
    }
    build_configuration config;
    config.k = 3;
    config.m = 1;
    config.verbose = false;
    config.tmp_dirname = directory.string();
    dictionary<aa_uint_kmer_t<uint64_t>, Offsets> dict;
    dict.build(input.string(), config);
    for (uint64_t i = 0; i < words.size(); ++i) {
        if (dict.lookup(words[i].data()).kmer_id != i) {
            throw std::runtime_error("incorrect heavy-bucket lookup for " + words[i]);
        }
    }
}


// BASIC graphs may retain both strands. A heavy bucket must still resolve the
// requested strand exactly, including when only its reverse complement exists.
template <typename Kmer, typename Offsets>
void check_forward_heavy(std::filesystem::path const& directory, uint64_t k) {
    auto reverse = [](std::string s) {
        std::reverse(s.begin(), s.end());
        for (char& c : s) c = c == 'A' ? 'T' : c == 'T' ? 'A' : c == 'C' ? 'G' : 'C';
        return s;
    };
    const std::string anchor = "TGACGATAAACCAACACGC";
    std::mt19937_64 random(20261104);
    auto flank = [&] {
        std::string s;
        for (uint64_t i = 0; i < (k - anchor.size()) / 2; ++i) s += "ACGT"[random() % 4];
        return s;
    };
    build_configuration config;
    config.k = k;
    config.m = anchor.size();
    config.allow_reverse_complement_duplicates = true;
    config.verbose = false;
    config.tmp_dirname = directory.string();
    hasher_type hasher(config.seed);
    auto anchor_kmer = util::string_to_uint_kmer<Kmer>(anchor.data(), anchor.size());
    uint64_t minimizer = util::canonical_mmer<Kmer>(uint64_t(anchor_kmer), anchor.size());
    std::vector<std::string> words;
    std::unordered_map<std::string, uint64_t> ids;
    for (uint64_t family = 0, attempts = 0; family < 512; ++attempts) {
        if (attempts > 10000) throw std::runtime_error("could not populate heavy bucket");
        auto word = flank() + anchor + flank();
        auto reversed = reverse(word);
        auto value = util::string_to_uint_kmer<Kmer>(word.data(), k);
        auto value_rc = value;
        value_rc.reverse_complement_inplace(k);
        if (ids.count(word) || ids.count(reversed)
            || util::compute_minimizer(value, value_rc, k, config.m, hasher).minimizer != minimizer)
            continue;
        ids.emplace(word, words.size());
        words.push_back(word);
        if (family % 4 < 2) {
            ids.emplace(reversed, words.size());
            words.push_back(reversed);
        }
        ++family;
    }
    auto input = directory / "forward.fa";
    {
        std::ofstream out(input);
        for (auto const& word : words) out << ">test\n" << word << '\n';
    }
    using dict_type = dictionary<Kmer, Offsets>;
    dict_type dict;
    dict.build(input.string(), config);
    uint64_t candidates = 0;
    dict.call_minimizer_positions(minimizer, [&](uint64_t) { ++candidates; });
    if (candidates <= (1ULL << constants::min_l))
        throw std::runtime_error("test did not produce a heavy bucket");
    auto verify = [&](dict_type const& index, std::string const& word) {
        for (bool with_rc : {false, true}) {
            auto result = index.lookup(word.data(), with_rc);
            bool present = ids.count(word) || (with_rc && ids.count(reverse(word)));
            if ((result.kmer_id != constants::invalid_uint64) != present)
                throw std::runtime_error("incorrect forward-only membership");
            if (!present) continue;
            std::string found(k, ' ');
            index.access(result.kmer_id, found.data());
            if (result.kmer_orientation == constants::backward_orientation) {
                if (!with_rc) throw std::runtime_error("forward-only lookup returned reverse strand");
                found = reverse(found);
            }
            if (found != word || (!with_rc && result.kmer_id != ids.at(word)))
                throw std::runtime_error("incorrect forward-only node or orientation");
        }
    };
    auto path = (directory / "forward.sshash").string();
    essentials::save(dict, path.c_str());
    dict_type loaded, mapped;
    essentials::load(loaded, path.c_str());
    essentials::mmap(mapped, path.c_str());
    for (auto word : words) {
        for (auto const* index : {&dict, &loaded, &mapped}) {
            verify(*index, word);
            verify(*index, reverse(word));
        }
        word[k / 2] = word[k / 2] == 'A' ? 'C' : 'A';
        verify(dict, word);
        verify(mapped, word);
    }
}

int main() {
    try {
        test_directory directory;
        // One record per distinct protein 3-mer produces dense minimizer buckets.
        std::vector<std::string> words;
        const std::string alphabet = "ACDEFGHIKLMNPQRSTVWY";
        for (char a : alphabet) for (char b : alphabet) for (char c : alphabet) {
            words.push_back(std::string{a, b, c});
        }
        check<decoded_offsets>(directory.path, words);
        check<encoded_offsets>(directory.path, words);
        check_forward_heavy<dna_uint_kmer_t<uint64_t>, decoded_offsets>(directory.path, 31);
        check_forward_heavy<dna_uint_kmer_t<uint64_t>, encoded_offsets>(directory.path, 31);
        check_forward_heavy<dna_uint_kmer_t<__uint128_t>, decoded_offsets>(directory.path, 41);
        check_forward_heavy<dna_uint_kmer_t<__uint128_t>, encoded_offsets>(directory.path, 41);
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
