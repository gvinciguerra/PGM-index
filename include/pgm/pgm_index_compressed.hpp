// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2019 Giorgio Vinciguerra.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include "sdsl.hpp"
#include "pgm_index.hpp"

namespace pgm {

/**
 * A space-efficient and compressed index that enables fast search operations on a sorted sequence of numbers.
 *
 * This is a variant of the @ref PGMIndex that internally uses compression to reduce the space of the index.
 *
 * @tparam K the type of the indexed keys
 * @tparam Epsilon controls the size of the returned search range
 * @tparam EpsilonRecursive controls the size of the search range in the internal structure
 * @tparam Floating the floating-point type to use for slopes
 */
template<typename K, size_t Epsilon, size_t EpsilonRecursive = 4, typename Floating = double>
class CompressedPGMIndex {
    static_assert(Epsilon > 0);
    struct CompressedLevel;

    size_t n;                             ///< The number of elements in the indexed data.
    K first_key;                          ///< The smallest element in the data.
    Floating root_slope;                  ///< The slope of the root segment.
    int64_t root_intercept;               ///< The intercept of the root segment.
    std::vector<Floating> slopes_table;   ///< The vector containing the slopes used by the segments in the index.
    std::vector<CompressedLevel> levels;  ///< The levels composing the compressed index.

    using floating_pair = std::pair<Floating, Floating>;
    using canonical_segment = typename internal::OptimalPiecewiseLinearModel<K, size_t>::CanonicalSegment;

public:

    static constexpr size_t epsilon_value = Epsilon;

    /**
     * Constructs an empty index.
     */
    CompressedPGMIndex() = default;

    /**
     * Constructs the compressed index on the given sorted vector.
     * @param data the vector of elements to be indexed, must be sorted
     */
    explicit CompressedPGMIndex(const std::vector<K> &data) : CompressedPGMIndex(data.begin(), data.end()) {}

    /**
     * Constructs the compressed index on the sorted elements in the range [first, last).
     * @param first, last the range containing the sorted elements to be indexed
     */
    template<typename Iterator>
    CompressedPGMIndex(Iterator first, Iterator last) : n(std::distance(first, last)) {
        if (n == 0)
            return;

        std::vector<size_t> levels_offsets({0});
        std::vector<canonical_segment> segments;
        segments.reserve(n / (Epsilon * Epsilon));

        auto ignore_last = *std::prev(last) == std::numeric_limits<K>::max(); // max is reserved for padding
        auto last_n = n - ignore_last;
        last -= ignore_last;

        // Build first level
        auto in_fun = [this, first](auto i) {
            auto x = first[i];
            if (i > 0 && i + 1u < n && x == first[i - 1] && x != first[i + 1] && x + 1 != first[i + 1])
                return std::pair<K, size_t>(x + 1, i);
            return std::pair<K, size_t>(x, i);
        };
        auto out_fun = [&, this](auto cs) { segments.emplace_back(cs); };
        last_n = internal::make_segmentation_par(last_n, Epsilon, in_fun, out_fun);
        levels_offsets.push_back(levels_offsets.back() + last_n);

        // Build upper levels
        while (EpsilonRecursive && last_n > 1) {
            auto offset = levels_offsets[levels_offsets.size() - 2];
            auto in_fun_rec = [&, this](auto i) { return std::pair<K, size_t>(segments[offset + i].get_first_x(), i); };
            last_n = internal::make_segmentation(last_n, EpsilonRecursive, in_fun_rec, out_fun);
            levels_offsets.push_back(levels_offsets.back() + last_n);
        }

        // Compress the slopes
        auto[tmp_table, map, intercepts] = merge_slopes(segments);
        slopes_table = tmp_table;

        // Build levels
        first_key = *first;
        if constexpr (EpsilonRecursive > 0) {
            auto root = *std::prev(levels_offsets.end(), 2);
            std::tie(root_slope, root_intercept) = segments[root].get_floating_point_segment(first_key);
        }

        levels.reserve(levels_offsets.size() - 2);
        for (auto i = EpsilonRecursive == 0 ? 1 : int(levels_offsets.size()) - 2; i > 0; --i) {
            auto l = levels_offsets[i - 1];
            auto r = levels_offsets[i];
            auto prev_level_size = i == 1 ? n : l - levels_offsets[i - 2];
            levels.emplace_back(segments.begin() + l, segments.begin() + r,
                                intercepts.begin() + l, intercepts.begin() + r,
                                map.begin() + l, map.begin() + r,
                                slopes_table, prev_level_size, *std::prev(last));
        }
    }

    /**
     * Returns the size of the index in bytes.
     * @return the size of the index in bytes
     */
    size_t size_in_bytes() const {
        size_t accum = 0;
        for (auto &l : levels)
            accum += l.size_in_bytes();
        return accum + slopes_table.size() * sizeof(Floating);
    }

    /**
     * Returns the approximate position and the range where @p key can be found.
     * @param key the value of the element to search for
     * @return a struct with the approximate position and bounds of the range
     */
    ApproxPos search(const K &key) const {
        auto k = std::max(first_key, key);

        if constexpr (EpsilonRecursive == 0) {
            auto &level = levels.front();
            auto it = std::upper_bound(level.keys.begin(), level.keys.begin() + level.size(), key);
            auto i = std::distance(level.keys.begin(), it);
            i = i == 0 ? 0 : i - 1;
            auto pos = std::min<size_t>(level(slopes_table, i, k), level.get_intercept(i + 1));
            auto lo = PGM_SUB_EPS(pos, Epsilon);
            auto hi = PGM_ADD_EPS(pos, Epsilon, n);
            return {pos, lo, hi};
        }

        auto p = int64_t(root_slope * (k - first_key)) + root_intercept;
        auto pos = std::min<size_t>(p > 0 ? size_t(p) : 0ull, levels.front().size());

        for (auto &level : levels) {
            auto lo = level.keys.begin() + PGM_SUB_EPS(pos, EpsilonRecursive + 1);

            static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(K);
            if constexpr (EpsilonRecursive <= linear_search_threshold) {
                for (; *std::next(lo) <= key; ++lo)
                    continue;
            } else {
                auto hi = level.keys.begin() + PGM_ADD_EPS(pos, EpsilonRecursive, level.size());
                auto it = std::upper_bound(lo, hi, k);
                lo == level.keys.begin() ? it : std::prev(it);
            }

            auto i = std::distance(level.keys.begin(), lo);
            pos = std::min<size_t>(level(slopes_table, i, k), level.get_intercept(i + 1));
        }

        auto lo = PGM_SUB_EPS(pos, Epsilon);
        auto hi = PGM_ADD_EPS(pos, Epsilon, n);
        return {pos, lo, hi};
    }

    /**
     * Returns the number of segments in the last level of the index.
     * @return the number of segments
     */
    size_t segments_count() const {
        return levels.back().size();
    }

    /**
     * Returns the number of levels of the index.
     * @return the number of levels of the index
     */
    size_t height() const {
        return levels.size() + 1;
    }

private:

    template<typename T, typename Cmp>
    static std::vector<size_t> sort_indexes(const std::vector<T> &v, Cmp cmp) {
        std::vector<size_t> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), cmp);
        return idx;
    }

    static std::tuple<std::vector<Floating>, std::vector<uint32_t>, std::vector<int64_t>>
    merge_slopes(const std::vector<canonical_segment> &segments) {
        std::vector<K> keys;
        std::vector<Floating> slopes_table;
        std::vector<int64_t> intercepts;
        intercepts.reserve(segments.size());
        slopes_table.reserve(segments.size());

        auto cmp = [&](auto i1, auto i2) { return segments[i1].get_slope_range() < segments[i2].get_slope_range(); };
        auto sorted_indexes = sort_indexes(segments, cmp);
        auto[current_min, current_max] = segments[sorted_indexes[0]].get_slope_range();

        std::vector<uint32_t> mapping(segments.size());
        mapping[sorted_indexes[0]] = 0;

        for (size_t i = 1; i < sorted_indexes.size(); ++i) {
            auto[min, max] = segments[sorted_indexes[i]].get_slope_range();
            if (min > current_max) {
                auto slope = 0.5 * (current_min + current_max);
                slopes_table.push_back(slope);
                current_min = min;
                current_max = max;
            } else {
                if (min > current_min)
                    current_min = min;
                if (max < current_max)
                    current_max = max;
            }
            mapping[sorted_indexes[i]] = uint32_t(slopes_table.size());
        }

        slopes_table.push_back(0.5 * (current_min + current_max));

        // Compute intercepts
        intercepts.reserve(segments.size());
        for (size_t i = 0; i < segments.size(); ++i) {
            auto[i_x, i_y] = segments[i].get_intersection();
            auto slope = slopes_table[mapping[i]];
            auto intercept = (int64_t) std::round(i_y - (i_x - segments[i].get_first_x()) * slope);
            intercepts.push_back(intercept);
        }

        return {slopes_table, mapping, intercepts};
    }
};

template<typename K, size_t Epsilon, size_t EpsilonRecursive, typename Floating>
struct CompressedPGMIndex<K, Epsilon, EpsilonRecursive, Floating>::CompressedLevel {
    std::vector<K> keys;                       ///< The keys of the segment in this level.
    sdsl::int_vector<> slopes_map;             ///< The ith element is an index into slopes_table.
    int64_t intercept_offset;                  ///< An offset to make the intercepts start from 0 in the bitvector.
    sdsl::sd_vector<> compressed_intercepts;   ///< The compressed bitvector storing the intercepts.
    sdsl::sd_vector<>::select_1_type sel1;     ///< The select1 succinct data structure on compressed_intercepts.

    template<typename IterK, typename IterI, typename IterM>
    CompressedLevel(IterK first_segment, IterK last_segment,
                    IterI first_intercept, IterI last_intercept,
                    IterM first_slope, IterM last_slope,
                    const std::vector<Floating> &slopes_table,
                    size_t prev_level_size,
                    K last_key)
        : keys(),
          intercept_offset(*first_intercept) {
        // If true, we need an extra segment to ensure that keys > *(last-1) are approximated to a position == n
        auto need_extra_segment = slopes_table[*std::prev(last_slope)] == 0;

        // Store keys
        keys.reserve(std::distance(first_segment, last_segment) + need_extra_segment + 1);
        for (auto it = first_segment; it != last_segment; ++it)
            keys.emplace_back(it->get_first_x());
        if (need_extra_segment)
            keys.emplace_back(last_key + 1);
        keys.emplace_back(std::numeric_limits<K>::max());

        // Compress and store intercepts
        sdsl::bit_vector intercept_bv(prev_level_size - intercept_offset + 1);
        intercept_bv[prev_level_size - intercept_offset] = true;
        for (auto it = first_intercept; it != last_intercept; ++it) {
            auto idx = std::min<int64_t>(prev_level_size - 1, *it) - intercept_offset;
            intercept_bv[idx] = true;
        }
        if (need_extra_segment)
            intercept_bv[prev_level_size + 1 - intercept_offset] = true;
        compressed_intercepts = sdsl::sd_vector<>(intercept_bv);
        sdsl::util::init_support(sel1, &compressed_intercepts);

        // Compress and store slopes_map
        size_t i = 0;
        size_t map_size = std::distance(first_slope, last_slope) + need_extra_segment;
        slopes_map = sdsl::int_vector<>(map_size, 0, sdsl::bits::hi(slopes_table.size() - 1) + 1);
        for (auto it = first_slope; it != last_slope; ++it)
            slopes_map[i++] = *it;
        if (need_extra_segment)
            slopes_map[slopes_map.size() - 1] = slopes_table[*std::prev(last_slope)];
    }

    inline size_t operator()(const std::vector<Floating> &slopes, size_t i, K k) const {
        auto pos = int64_t(get_slope(slopes, i) * (k - keys[i])) + get_intercept(i);
        return pos > 0 ? size_t(pos) : 0ull;
    }

    inline Floating get_slope(const std::vector<Floating> &slopes, size_t i) const {
        return slopes[slopes_map[i]];
    }

    inline int64_t get_intercept(size_t i) const {
        return intercept_offset + int64_t(sel1(i + 1));
    }

    inline size_t size() const {
        return keys.size() - 1;
    }

    inline size_t size_in_bytes() const {
        return keys.size() * sizeof(K) + slopes_map.bit_size() / 8 + sdsl::size_in_bytes(compressed_intercepts);
    }
};

}