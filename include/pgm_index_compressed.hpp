// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2019 Giorgio Vinciguerra.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#pragma once

#include <sdsl/bit_vectors.hpp>
#include "pgm_index.hpp"

/**
 * A space-efficient and compressed index that finds the position of a sought key within a radius of @p Error.
 *
 * This is a variant of the @ref PGMIndex that internally uses compression to reduce the space of the index.
 *
 * @tparam K the type of the indexed elements
 * @tparam Error the maximum error allowed in the last level of the index
 * @tparam RecursiveError the maximum error allowed in the upper levels of the index
 * @tparam Floating the floating-point type to use for slopes
 */
template<typename K, size_t Error, size_t RecursiveError = 8, typename Floating = double>
class CompressedPGMIndex {
    struct CompressedLevel;

    size_t n;                             ///< The number of elements in the indexed data.
    K root_key;                           ///< The key of the root segment.
    Floating root_slope;                  ///< The slope of the root segment.
    Floating root_intercept;              ///< The intercept of the root segment.
    K first_key;                          ///< The smallest element in the data.
    K last_key;                           ///< The largest element in the data.
    std::vector<Floating> slopes_table;   ///< The vector containing the slopes used by the segments in the index.
    std::vector<CompressedLevel> levels;  ///< The levels composing the compressed index.

    struct CompressedLevel {
        std::vector<K> keys;                       ///< The keys of the segment in this level.
        sdsl::int_vector<> slopes_map;             ///< The ith element is an index into slopes_table.
        std::vector<Floating> &slopes_table;       ///< A reference to the vector containing the slopes.
        int64_t intercept_offset;                  ///< An offset to make the intercepts start from 0 in the bitvector.
        sdsl::sd_vector<> compressed_intercepts;   ///< The compressed bitvector storing the intercepts.
        sdsl::sd_vector<>::select_1_type sel1;     ///< The select1 succinct data structure on compressed_intercepts.

        template<typename IterK, typename IterI, typename IterM>
        CompressedLevel(IterK first_key, IterK last_key,
                        IterI first_intercept, IterI last_intercept,
                        IterM first_slope, IterM last_slope,
                        std::vector<Floating> &slopes_table)
            : keys(first_key, last_key), slopes_table(slopes_table), intercept_offset(*first_intercept) {
            sdsl::bit_vector intercept_bv(*(last_intercept - 1) - intercept_offset + 1);
            for (auto it = first_intercept; it != last_intercept; ++it)
                intercept_bv[*it - intercept_offset] = 1;
            compressed_intercepts = sdsl::sd_vector<>(intercept_bv);
            sdsl::util::init_support(sel1, &compressed_intercepts);

            size_t i = 0;
            size_t map_size = std::distance(first_slope, last_slope);
            slopes_map = sdsl::int_vector<>(map_size, 0, sdsl::bits::hi(slopes_table.size() - 1) + 1);
            for (auto it = first_slope; it != last_slope; ++it)
                slopes_map[i++] = *it;
        }

        inline size_t operator()(size_t i, K k) const {
            Floating pos = get_slope(i) * (k - keys[i]) + get_intercept(i);
            return pos > Floating(0) ? pos : 0ul;
        }

        inline Floating get_slope(size_t i) const {
            return slopes_table[slopes_map[i]];
        }

        inline Floating get_intercept(size_t i) const {
            return intercept_offset + int64_t(sel1(i + 1));
        }

        inline size_t size() const {
            return keys.size();
        }

        inline size_t size_in_bytes() const {
            return keys.size() * sizeof(K) + slopes_map.bit_size() / 8 + sdsl::size_in_bytes(compressed_intercepts);
        }
    };

public:

    /**
     * Constructs the compressed index on the given sorted data.
     * @param data the vector of keys, must be sorted
     */
    explicit CompressedPGMIndex(std::vector<K> &data) : CompressedPGMIndex(data.begin(), data.end()) {}

    /**
     * Constructs the compressed index on the sorted data in the range [first, last).
     * @param first, last the range containing the sorted elements to be indexed
     */
    template<typename Iterator>
    CompressedPGMIndex(Iterator first, Iterator last) : n(std::distance(first, last)) {
        assert(std::is_sorted(first, last));
        if (n == 0)
            return;

        first_key = *first;
        last_key = *std::prev(last);

        std::vector<size_t> levels_boundaries({0});
        std::vector<K> tmp_keys;
        std::vector<std::pair<Floating, Floating>> tmp_ranges;
        std::vector<std::pair<Floating, Floating>> tmp_intersections;
        tmp_keys.reserve(2048);
        tmp_ranges.reserve(2048);
        tmp_intersections.reserve(2048);
        bool needs_more_levels = true;

        // Iterative construction of the levels
        while (needs_more_levels) {
            bool first_level = levels_boundaries.size() == 1;
            size_t error = first_level ? Error : RecursiveError;
            size_t offset = first_level ? 0 : levels_boundaries[levels_boundaries.size() - 2];
            size_t last_n = first_level ? n : tmp_keys.size();
            auto it = first_level ? first : tmp_keys.begin() + offset;
            auto key = *it;

            // Compute segmentation for this level
            OptimalPiecewiseLinearModel<K, size_t> algorithm(error);
            algorithm.add_point(*it, 0);
            ++it;

            for (size_t i = offset + 1; i < last_n; ++i, ++it) {
                if (*it == *std::prev(it))
                    continue;

                if (!algorithm.add_point(*it, i - offset)) {
                    auto cs = algorithm.get_segment();
                    tmp_keys.emplace_back(key);
                    tmp_ranges.emplace_back(cs.get_slope_range());
                    tmp_intersections.emplace_back(cs.get_intersection());
                    key = *it;
                    --i;
                    --it;
                }
            }

            // Last segment
            auto cs = algorithm.get_segment();
            tmp_keys.emplace_back(key);
            tmp_ranges.emplace_back(cs.get_slope_range());
            tmp_intersections.emplace_back(cs.get_intersection());

            needs_more_levels = tmp_keys.size() - levels_boundaries.back() > 1;
            levels_boundaries.push_back(tmp_keys.size());
        }

        // Compress the slopes
        auto[tmp_table, tmp_map, tmp_intercepts] = merge_slopes(tmp_ranges, tmp_intersections, tmp_keys);
        slopes_table = tmp_table;

        // Build levels
        auto[i_x, i_y] = tmp_intersections.back();
        root_key = tmp_keys.back();
        root_slope = 0.5 * (tmp_ranges.back().first + tmp_ranges.back().second);
        root_intercept = std::floor(i_y - (i_x - root_key) * root_slope);

        levels.reserve(levels_boundaries.size() - 2);
        for (int i = levels_boundaries.size() - 2; i > 0; --i) {
            auto l = levels_boundaries[i - 1];
            auto r = levels_boundaries[i];
            levels.emplace_back(tmp_keys.begin() + l, tmp_keys.begin() + r,
                                tmp_intercepts.begin() + l, tmp_intercepts.begin() + r,
                                tmp_map.begin() + l, tmp_map.begin() + r, slopes_table);
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
     * Returns the approximate position of a key.
     * @param key the value of the element to search for
     * @return a struct with the approximate position
     */
    ApproxPos find_approximate_position(K key) const {
        if (key <= first_key)
            return {0, 0, 1};
        if (key >= last_key)
            return {n - 1, n - 1, n};

        size_t approx_pos = std::max<Floating>(0, root_intercept + root_slope * (key - root_key));
        size_t pos = 0;

        for (auto &level : levels) {
            auto lo = level.keys.begin() + SUB_ERR(approx_pos, RecursiveError + 1);
            auto hi = level.keys.begin() + ADD_ERR(approx_pos, RecursiveError + 1, level.size());

            static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(K);
            if constexpr (RecursiveError <= linear_search_threshold) {
                for (; lo < hi && *lo <= key; ++lo);
                --lo;
            } else {
                auto it = std::upper_bound(lo, hi, key);
                lo == level.keys.begin() ? it : std::prev(it);
            }
            pos = std::distance(level.keys.begin(), lo);

            approx_pos = level(pos, key);
            if (pos + 1 < level.size())
                approx_pos = std::min(approx_pos, (size_t) level.get_intercept(pos + 1));

            assert(level.keys[pos] <= key);
            assert(pos + 1 >= level.size() || level.keys[pos + 1] > key);
        }

        auto p = levels.back()(pos, key);
        auto lo = SUB_ERR(p, Error);
        auto hi = ADD_ERR(p, Error + 1, n);
        p = p >= n ? n - 1 : p;

        return {p, lo, hi};
    }

    /**
     * Returns the number of segments in the last level of the index.
     * @return the number of segments
     */
    size_t segments_count() const {
        return levels.back().size();
    }

    /**
     * Returns the number of levels in the index.
     * @return the number of levels in the index
     */
    size_t height() const {
        return levels.size() + 1;
    }

private:

    template<typename T>
    static std::vector<size_t> sort_indexes(const std::vector<T> &v) {
        std::vector<size_t> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
        return idx;
    }

    static std::tuple<std::vector<Floating>, std::vector<uint32_t>, std::vector<int64_t>>
    merge_slopes(const std::vector<std::pair<Floating, Floating>> &slope_ranges,
                 const std::vector<std::pair<Floating, Floating>> &intersections,
                 const std::vector<K> &keys) {
        std::vector<Floating> slopes_table;
        std::vector<int64_t> intercepts;
        intercepts.reserve(intersections.size());
        slopes_table.reserve(slope_ranges.size());
        auto sorted_indexes = sort_indexes(slope_ranges);
        auto[current_min, current_max] = slope_ranges[sorted_indexes[0]];
        std::vector<uint32_t> tmp_mapping(slope_ranges.size());
        tmp_mapping[sorted_indexes[0]] = 0;

        for (size_t i = 1; i < sorted_indexes.size(); ++i) {
            auto[min, max] = slope_ranges[sorted_indexes[i]];
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
            tmp_mapping[sorted_indexes[i]] = uint32_t(slopes_table.size());
        }

        slopes_table.push_back(0.5 * (current_min + current_max));

        // Compute intercepts
        intercepts.reserve(intersections.size());
        for (size_t i = 0; i < intersections.size(); ++i) {
            auto[i_x, i_y] = intersections[i];
            auto slope = slopes_table[tmp_mapping[i]];
            auto intercept = std::floor(i_y - (i_x - keys[i]) * slope);
            intercepts.push_back(int64_t(intercept));
        }

        return {slopes_table, tmp_mapping, intercepts};
    }

};