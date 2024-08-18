#include "dictionary.hpp"
#include "skew_index.hpp"

namespace sshash {

uint64_t skew_index::print_info() const {
    uint64_t num_partitions = mphfs.size();
    uint64_t lower = 1ULL << min_log2;
    uint64_t upper = 2 * lower;
    uint64_t num_kmers_in_skew_index = 0;
    for (uint64_t partition_id = 0; partition_id != num_partitions; ++partition_id) {
        uint64_t n = mphfs[partition_id].num_keys();
        assert(n == positions[partition_id].size());
        std::cout << "num_kmers belonging to buckets of size > " << lower << " and <= " << upper
                  << ": " << n << "; ";
        std::cout << "bits/kmer = " << static_cast<double>(mphfs[partition_id].num_bits()) / n
                  << " (mphf) + " << (positions[partition_id].bytes() * 8.0) / n
                  << " (positions)\n";
        num_kmers_in_skew_index += n;
        lower = upper;
        upper = 2 * lower;
    }
    return num_kmers_in_skew_index;
}

double bits_per_kmer_formula(uint64_t k, /* kmer length */
                             uint64_t m, /* minimizer length */
                             uint64_t n, /* num. kmers */
                             uint64_t M) /* num. strings in SPSS */
{
    /*
      Caveats:
      1. we assume an alphabet of size 4
      2. this assumes a random minimizer scheme, so num. super-kmers is ~ 2n/(k-m+2)
      3. we neglect lower order terms and skew index space
      4. no canonical parsing
    */

    assert(k > 0);
    assert(k >= m);

    const uint64_t N = n + M * (k - 1);  // num. characters in SPSS

    /* summing (M-1) provides an upper bound to the num. of super-kmers */
    double Z = (2.0 * n) / (k - m + 2) + (M - 1);

    /* A cache line is 64 B = 512 bits -->
       max window_size that fits in a cache line is 512/2 = 256
       assuming a 2-bit encoded stream. */
    const uint64_t window_size = 1; /* 256; */

    double num_bits =
        2 * N + Z * (5.0 + std::ceil(std::log2(std::ceil(static_cast<double>(N) / window_size)))) +
        M * (2.0 + std::ceil(std::log2(static_cast<double>(N) / M)));

    return num_bits / n;
}

double perc(uint64_t amount, uint64_t total) { return (amount * 100.0) / total; }

void dictionary::print_space_breakdown() const {
    const uint64_t num_bytes = (num_bits() + 7) / 8;
    std::cout << "total index size: " << num_bytes << " [B] -- "
              << essentials::convert(num_bytes, essentials::MB) << " [MB]" << '\n';
    std::cout << "SPACE BREAKDOWN:\n";
    std::cout << "  minimizers: " << static_cast<double>(m_minimizers.num_bits()) / size()
              << " [bits/kmer] ("
              << static_cast<double>(m_minimizers.num_bits()) / m_minimizers.size()
              << " [bits/key]) -- " << perc(m_minimizers.num_bits(), num_bits()) << "%\n";
    std::cout << "  pieces: " << static_cast<double>(m_buckets.pieces.num_bits()) / size()
              << " [bits/kmer] -- " << perc(m_buckets.pieces.num_bits(), num_bits()) << "%\n";
    std::cout << "  num_super_kmers_before_bucket: "
              << static_cast<double>(m_buckets.num_super_kmers_before_bucket.num_bits()) / size()
              << " [bits/kmer] -- "
              << perc(m_buckets.num_super_kmers_before_bucket.num_bits(), num_bits()) << "%\n";
    std::cout << "  offsets: " << static_cast<double>(8 * m_buckets.offsets.bytes()) / size()
              << " [bits/kmer] -- " << perc(8 * m_buckets.offsets.bytes(), num_bits()) << "%\n";
    std::cout << "  strings: " << static_cast<double>(8 * m_buckets.strings.bytes()) / size()
              << " [bits/kmer] -- " << perc(8 * m_buckets.strings.bytes(), num_bits()) << "%\n";
    std::cout << "  skew_index: " << static_cast<double>(m_skew_index.num_bits()) / size()
              << " [bits/kmer] -- " << perc(m_skew_index.num_bits(), num_bits()) << "%\n";
    std::cout << "  weights: " << static_cast<double>(m_weights.num_bits()) / size()
              << " [bits/kmer] -- " << perc(m_weights.num_bits(), num_bits()) << "%\n";
    m_weights.print_space_breakdown(size());
    std::cout << "  --------------\n";
    std::cout << "  total: " << static_cast<double>(num_bits()) / size() << " [bits/kmer]"
              << std::endl;

    std::cout << "  Close-form formula: " << bits_per_kmer_formula(k(), m(), size(), num_contigs())
              << " [bits/kmer]" << std::endl;
}

void dictionary::print_info() const {
    std::cout << "=== dictionary info:\n";
    std::cout << "num_kmers = " << size() << '\n';
    std::cout << "k = " << k() << '\n';
    std::cout << "num_minimizers = " << m_minimizers.size() << std::endl;
    std::cout << "m = " << m() << '\n';
    std::cout << "canonicalized = " << (canonicalized() ? "true" : "false") << '\n';
    std::cout << "weighted = " << (weighted() ? "true" : "false") << '\n';

    std::cout << "num_super_kmers = " << m_buckets.offsets.size() << '\n';
    std::cout << "num_pieces = " << m_buckets.pieces.size() << " (+"
              << (2.0 * m_buckets.pieces.size() * (k() - 1)) / size() << " [bits/kmer])" << '\n';
    std::cout << "bits_per_offset = ceil(log2(" << m_buckets.strings.size() / 2
              << ")) = " << std::ceil(std::log2(m_buckets.strings.size() / 2)) << '\n';
    uint64_t num_kmers_in_skew_index = m_skew_index.print_info();
    std::cout << "num_kmers_in_skew_index " << num_kmers_in_skew_index << "("
              << (num_kmers_in_skew_index * 100.0) / size() << "%)" << std::endl;

    print_space_breakdown();
}

}  // namespace sshash