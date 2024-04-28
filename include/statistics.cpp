#include "dictionary.hpp"
#include "buckets_statistics.hpp"
namespace sshash {

void dictionary::compute_statistics() const {
    uint64_t num_kmers = size();
    uint64_t num_minimizers = m_minimizers.size();
    uint64_t num_super_kmers = m_buckets.offsets.size();

    buckets_statistics buckets_stats(num_minimizers, num_kmers, num_super_kmers);

    std::cout << "computing buckets statistics..." << std::endl;

    for (uint64_t bucket_id = 0; bucket_id != num_minimizers; ++bucket_id) {
        auto [begin, end] = m_buckets.locate_bucket(bucket_id);
        uint64_t num_super_kmers_in_bucket = end - begin;
        buckets_stats.add_num_super_kmers_in_bucket(num_super_kmers_in_bucket);
        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
            uint64_t offset = m_buckets.offsets.access(super_kmer_id);
            auto [_, contig_end] = m_buckets.offset_to_id(offset, m_k);
            (void)_;
            bit_vector_iterator bv_it(m_buckets.strings, 2 * offset);
            uint64_t window_size = std::min<uint64_t>(m_k - m_m + 1, contig_end - offset - m_k + 1);
            uint64_t prev_minimizer = constants::invalid_uint64;
            uint64_t w = 0;
            for (; w != window_size; ++w) {
                uint64_t kmer = bv_it.read_and_advance_by_two(2 * m_k);
                auto [minimizer, pos] = util::compute_minimizer_pos(kmer, m_k, m_m, m_seed);
                if (m_canonical_parsing) {
                    uint64_t kmer_rc = util::compute_reverse_complement(kmer, m_k);
                    auto [minimizer_rc, pos_rc] =
                        util::compute_minimizer_pos(kmer_rc, m_k, m_m, m_seed);
                    if (minimizer_rc < minimizer) {
                        minimizer = minimizer_rc;
                        pos = pos_rc;
                    }
                }
                if (prev_minimizer != constants::invalid_uint64 and minimizer != prev_minimizer) {
                    break;
                }
                prev_minimizer = minimizer;
            }
            buckets_stats.add_num_kmers_in_super_kmer(num_super_kmers_in_bucket, w);
        }
    }
    buckets_stats.print_full();
    std::cout << "DONE" << std::endl;
}
inline bool equal(const std::vector<std::string>& input1, const std::vector<std::string>& input2) {
    if(input1.size() != input2.size()){
        return false;
    }
    for(size_t i = 0; i < input1.size(); i++){
        if(input1[i] != input2[i]){
            return false;
        }
    }
    return true; 
}

void uint_kmer_to_last_char(kmer_t x, char* str, uint64_t k){
    x >>= 2*(k-1);
    str[0] = util::uint64_to_char(x & 3);
}

std::vector<bool> dictionary::build_superkmer_bv(const std::function <bool (std::string_view)> &monochromatic_labels) const {
    //uint64_t num_kmers = size();
    uint64_t num_minimizers = m_minimizers.size();
    uint64_t num_super_kmers = m_buckets.offsets.size();

    std::cout << "building super kmer mask..." << std::endl;
    std::vector<bool> non_mono_superkmer (num_super_kmers, false);
    //size_t superkmer_idx = 0;
    uint64_t one_pc_buckets = std::ceil(num_minimizers/100.0);
    std::cout<<"iterating through buckets :\n";
    uint64_t num_threads = std::thread::hardware_concurrency();
    if (num_minimizers < num_threads) num_threads = num_minimizers;
    std::mutex vec_mutex;
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (uint64_t bucket_id = 0; bucket_id != num_minimizers; ++bucket_id) {
        if(bucket_id%one_pc_buckets==0)std::cout << bucket_id/one_pc_buckets <<"%" << '\r'<< std::flush;
	    auto [begin, end] = m_buckets.locate_bucket(bucket_id);
        //uint64_t num_super_kmers_in_bucket = end - begin;
        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
	    uint64_t offset = m_buckets.offsets.access(super_kmer_id);
            auto [_, contig_end] = m_buckets.offset_to_id(offset, m_k);
            (void)_;
            bit_vector_iterator bv_it(m_buckets.strings, 2 * offset);
            uint64_t window_size = std::min<uint64_t>(m_k - m_m + 1, contig_end - offset - m_k + 1);
            uint64_t prev_minimizer = constants::invalid_uint64;
            std::vector<std::string> prev_labels = {};
            std::string superkmer = "";
            size_t num_chars_to_add = m_k;
            for (uint64_t w = 0; w != window_size; ++w) {
                kmer_t kmer = bv_it.read_and_advance_by_two(2 * m_k);
                auto [minimizer, pos] = util::compute_minimizer_pos(kmer, m_k, m_m, m_seed);
                
                //check if superkmer has ended
                if (prev_minimizer != constants::invalid_uint64 and minimizer != prev_minimizer) {
                    break;
                }
                prev_minimizer = minimizer;

                //get kmers into superkmer
                std::string chars_to_add(num_chars_to_add, '_');
                if(num_chars_to_add > 1){
                    util::uint_kmer_to_string(kmer, &chars_to_add[0], m_k);
		    num_chars_to_add = 1;
                }else{
                    uint_kmer_to_last_char(kmer, &chars_to_add[0], m_k);
                }
                superkmer += chars_to_add;
                
            }   
            // get labels and compare them to previous ones
            if(!monochromatic_labels(superkmer)){
                vec_mutex.lock();
		non_mono_superkmer[super_kmer_id] = true;
		vec_mutex.unlock();
            }
            
        }
    }
    std::cout << "DONE" << std::endl;

    return non_mono_superkmer;
}

}  // namespace sshash
