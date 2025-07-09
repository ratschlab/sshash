#pragma once

#include "external/pthash/external/bits/include/elias_fano.hpp"

#include "util.hpp"
#include "kmer_iterator.hpp"

namespace sshash {

template <class kmer_t>
struct buckets  //
{
    lookup_result offset_to_id(const uint64_t offset, const uint64_t k) const {
        auto p = pieces.locate(offset);
        uint64_t contig_id = p.first.pos;
        uint64_t contig_begin = p.first.val;
        uint64_t contig_end = p.second.val;

        /* The following facts hold. */
        assert(offset >= contig_id * (k - 1));
        assert(contig_begin <= offset);
        assert(offset < contig_end);
        /****************************/

        uint64_t absolute_kmer_id = offset - contig_id * (k - 1);
        uint64_t relative_kmer_id = offset - contig_begin;
        uint64_t contig_length = contig_end - contig_begin;
        assert(contig_length >= k);
        uint64_t contig_size = contig_length - k + 1;

        lookup_result res;
        res.kmer_id = absolute_kmer_id;
        res.kmer_id_in_contig = relative_kmer_id;
        res.contig_id = contig_id;
        res.contig_size = contig_size;
        assert(contig_begin == res.contig_begin(k));
        assert(contig_end == res.contig_end(k));

        return res;
    }

    /* Return where the contig begins and ends in strings. */
    std::pair<uint64_t, uint64_t>  // [begin, end)
    contig_offsets(const uint64_t contig_id) const {
        uint64_t begin = pieces.access(contig_id);
        uint64_t end = pieces.access(contig_id + 1);
        assert(end > begin);
        return {begin, end};
    }

    kmer_t contig_prefix(const uint64_t contig_id, const uint64_t k) const {
        uint64_t contig_begin = pieces.access(contig_id);
        return util::read_kmer_at<kmer_t>(strings, k - 1, kmer_t::bits_per_char * contig_begin);
    }

    kmer_t contig_suffix(const uint64_t contig_id, const uint64_t k) const {
        uint64_t contig_end = pieces.access(contig_id + 1);
        return util::read_kmer_at<kmer_t>(strings, k - 1,
                                          kmer_t::bits_per_char * (contig_end - k + 1));
    }

    std::pair<uint64_t, uint64_t> locate_bucket(const uint64_t bucket_id) const {
        uint64_t begin = num_super_kmers_before_bucket.access(bucket_id) + bucket_id;
        uint64_t end = num_super_kmers_before_bucket.access(bucket_id + 1) + bucket_id + 1;
        assert(begin < end);
        return {begin, end};
    }

    lookup_result lookup(uint64_t bucket_id, kmer_t target_kmer,                               //
                         uint64_t target_minimizer,                                            //
                         const uint64_t k, const uint64_t m, hasher_type const& hasher) const  //
    {
        auto [begin, end] = locate_bucket(bucket_id);
        return lookup(begin, end, target_kmer, k, m, target_minimizer, hasher);
    }

    lookup_result lookup(uint64_t begin, uint64_t end, kmer_t target_kmer,                     //
                         uint64_t target_minimizer,                                            //
                         const uint64_t k, const uint64_t m, hasher_type const& hasher) const  //
    {
        { /* check minimizer first */
            uint64_t offset = offsets.access(begin);
            auto read_kmer = util::read_kmer_at<kmer_t>(strings, k, kmer_t::bits_per_char * offset);
            uint64_t minimizer = util::compute_minimizer(read_kmer, k, m, hasher);
            if (minimizer != target_minimizer) {
                auto res = lookup_result();
                res.minimizer_found = false;
                return res;
            }
        }

        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
            auto res = lookup_in_super_kmer(super_kmer_id, target_kmer, k, m);
            if (res.kmer_id != constants::invalid_uint64) {
                assert(is_valid(res));
                return res;
            }
        }

        return lookup_result();
    }

    lookup_result lookup_in_super_kmer(uint64_t super_kmer_id, kmer_t target_kmer,  //
                                       const uint64_t k, const uint64_t m) const    //
    {
        uint64_t offset = offsets.access(super_kmer_id);
        auto res = offset_to_id(offset, k);
        kmer_iterator<kmer_t> it(strings, k, kmer_t::bits_per_char * offset);
        uint64_t window_size = std::min<uint64_t>(k - m + 1, res.contig_end(k) - offset - k + 1);
        for (uint64_t w = 0; w != window_size; ++w) {
            auto read_kmer = it.get();
            if (read_kmer == target_kmer) {
                res.kmer_id += w;
                res.kmer_id_in_contig += w;
                assert(is_valid(res));
                return res;
            }
            it.next();
        }
        return lookup_result();
    }

    lookup_result lookup_canonical(uint64_t bucket_id, kmer_t target_kmer, kmer_t target_kmer_rc,
                                   uint64_t target_minimizer,           //
                                   const uint64_t k, const uint64_t m,  //
                                   hasher_type const& hasher) const     //
    {
        auto [begin, end] = locate_bucket(bucket_id);
        return lookup_canonical(begin, end, target_kmer, target_kmer_rc,  //
                                target_minimizer, k, m, hasher);
    }

    lookup_result lookup_canonical_in_super_kmer(uint64_t super_kmer_id,                     //
                                                 kmer_t target_kmer, kmer_t target_kmer_rc,  //
                                                 const uint64_t k, const uint64_t m) const   //
    {
        uint64_t offset = offsets.access(super_kmer_id);
        auto res = offset_to_id(offset, k);
        kmer_iterator<kmer_t> it(strings, k, kmer_t::bits_per_char * offset);
        uint64_t window_size = std::min<uint64_t>(k - m + 1, res.contig_end(k) - offset - k + 1);
        for (uint64_t w = 0; w != window_size; ++w) {
            auto read_kmer = it.get();
            if (read_kmer == target_kmer) {
                res.kmer_id += w;
                res.kmer_id_in_contig += w;
                res.kmer_orientation = constants::forward_orientation;
                assert(is_valid(res));
                return res;
            }
            if (read_kmer == target_kmer_rc) {
                res.kmer_id += w;
                res.kmer_id_in_contig += w;
                res.kmer_orientation = constants::backward_orientation;
                assert(is_valid(res));
                return res;
            }
            it.next();
        }
        return lookup_result();
    }

    lookup_result lookup_canonical(uint64_t begin, uint64_t end,               //
                                   kmer_t target_kmer, kmer_t target_kmer_rc,  //
                                   uint64_t target_minimizer,                  //
                                   const uint64_t k, const uint64_t m,         //
                                   hasher_type const& hasher) const            //
    {
        { /* check minimizer first */
            uint64_t offset = offsets.access(begin);
            auto read_kmer = util::read_kmer_at<kmer_t>(strings, k, kmer_t::bits_per_char * offset);
            kmer_t read_kmer_rc = read_kmer;
            read_kmer_rc.reverse_complement_inplace(k);
            uint64_t minimizer = std::min(util::compute_minimizer(read_kmer, k, m, hasher),
                                          util::compute_minimizer(read_kmer_rc, k, m, hasher));
            if (minimizer != target_minimizer) {
                auto res = lookup_result();
                res.minimizer_found = false;
                return res;
            }
        }

        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
            uint64_t offset = offsets.access(super_kmer_id);
            auto res = offset_to_id(offset, k);
            kmer_iterator<kmer_t> it(strings, k, kmer_t::bits_per_char * offset);
            uint64_t window_size =
                std::min<uint64_t>(k - m + 1, res.contig_end(k) - offset - k + 1);
            for (uint64_t w = 0; w != window_size; ++w) {
                auto read_kmer = it.get();
                if (read_kmer == target_kmer) {
                    res.kmer_id += w;
                    res.kmer_id_in_contig += w;
                    res.kmer_orientation = constants::forward_orientation;
                    assert(is_valid(res));
                    return res;
                }
                if (read_kmer == target_kmer_rc) {
                    res.kmer_id += w;
                    res.kmer_id_in_contig += w;
                    res.kmer_orientation = constants::backward_orientation;
                    assert(is_valid(res));
                    return res;
                }
                it.next();
            }
        }

        return lookup_result();
    }

    uint64_t id_to_offset(const uint64_t id, const uint64_t k) const {
        constexpr uint64_t linear_scan_threshold = 8;
        uint64_t lo = 0;
        uint64_t hi = pieces.size() - 1;
        assert(pieces.access(0) == 0);
        while (lo < hi) {
            if (hi - lo <= linear_scan_threshold) {
                for (; lo < hi; ++lo) {
                    uint64_t val = pieces.access(lo) - lo * (k - 1);
                    if (val > id) break;
                }
                break;
            }
            uint64_t mid = lo + (hi - lo) / 2;
            uint64_t val = pieces.access(mid);
            assert(val >= mid * (k - 1));
            if (id <= val - mid * (k - 1)) {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        if (lo < pieces.size() and pieces.access(lo) - lo * (k - 1) > id) --lo;
        return id + lo * (k - 1);
    }

    void access(const uint64_t kmer_id, char* string_kmer, const uint64_t k) const {
        uint64_t offset = id_to_offset(kmer_id, k);
        auto read_kmer = util::read_kmer_at<kmer_t>(strings, k, kmer_t::bits_per_char * offset);
        util::uint_kmer_to_string(read_kmer, string_kmer, k);
    }

    struct iterator {
        iterator() {}

        iterator(buckets const* ptr,                                        //
                 const uint64_t begin_kmer_id, const uint64_t end_kmer_id,  // [begin,end)
                 const uint64_t k)
            : m_buckets(ptr)
            , m_begin_kmer_id(begin_kmer_id)
            , m_end_kmer_id(end_kmer_id)
            , m_k(k)
            , m_it(ptr->strings, m_k)  //
        {
            m_offset = m_buckets->id_to_offset(m_begin_kmer_id, k);
            auto [pos, piece_end] = m_buckets->pieces.next_geq(m_offset);
            if (piece_end == m_offset) pos += 1;
            m_pieces_it = m_buckets->pieces.get_iterator_at(pos);
            next_piece();
            m_ret.second.resize(m_k, 0);
        }

        bool has_next() const { return m_begin_kmer_id != m_end_kmer_id; }

        std::pair<uint64_t, std::string> next() {
            if (m_offset == m_next_offset - m_k + 1) {
                m_offset = m_next_offset;
                next_piece();
            }
            m_ret.first = m_begin_kmer_id;
            if (m_clear) {
                util::uint_kmer_to_string(m_it.get(), m_ret.second.data(), m_k);
                assert(kmer_t::bits_per_char * m_offset == m_it.position());
                m_it.at(kmer_t::bits_per_char * (m_offset + m_k));
            } else {
                memmove(m_ret.second.data(), m_ret.second.data() + 1, m_k - 1);
                m_ret.second[m_k - 1] = kmer_t::uint64_to_char(m_it.get_next_char());
            }
            m_clear = false;
            ++m_begin_kmer_id;
            ++m_offset;
            return m_ret;
        }

    private:
        std::pair<uint64_t, std::string> m_ret;
        buckets const* m_buckets;
        uint64_t m_begin_kmer_id, m_end_kmer_id;
        uint64_t m_k;
        uint64_t m_offset;
        uint64_t m_next_offset;
        kmer_iterator<kmer_t> m_it;
        bits::elias_fano<true, false>::iterator m_pieces_it;
        bool m_clear;

        void next_piece() {
            m_it.at(kmer_t::bits_per_char * m_offset);
            m_next_offset = m_pieces_it.value();
            assert(m_next_offset > m_offset);
            m_clear = true;
            m_pieces_it.next();
        }
    };

    iterator at(const uint64_t begin_kmer_id, const uint64_t end_kmer_id, const uint64_t k) const {
        return iterator(this, begin_kmer_id, end_kmer_id, k);
    }

    uint64_t num_bits() const {
        return 8 * (pieces.num_bytes() + num_super_kmers_before_bucket.num_bytes() +
                    offsets.num_bytes() + strings.num_bytes());
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

    bits::elias_fano<true, false> pieces;
    bits::elias_fano<false, false> num_super_kmers_before_bucket;
    bits::compact_vector offsets;
    bits::bit_vector strings;

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.pieces);
        visitor.visit(t.num_super_kmers_before_bucket);
        visitor.visit(t.offsets);
        visitor.visit(t.strings);
    }

    bool is_valid(lookup_result res) const {
        return res.contig_size != constants::invalid_uint64 and  //
               res.kmer_id_in_contig < res.contig_size and       //
               res.contig_id < pieces.size();                    //
    }
};

}  // namespace sshash