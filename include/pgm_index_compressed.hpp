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
#include "packed_vector.hpp"

/**
 * A space-efficient and compressed index that finds the position of a sought key within a radius of @p Error.
 * @tparam K the type of the indexed elements
 * @tparam Error the maximum allowed error in the last level of the index
 * @tparam RecursiveError the maximum allowed error in the upper levels of the index
 * @tparam Floating the floating-point type to use for slopes
 * @tparam CompressedBV the type of the compressed bitvector storing the intercepts
 */
template<typename K, size_t Error, size_t RecursiveError = 16, typename Floating = double,
    typename CompressedBV = sdsl::sd_vector<>>
class CompressedPGMIndex {
    size_t data_size;                   ///< The number of elements in the indexed data.
    std::vector<Floating> slopes_table; ///< The vector containing the slopes used by the segments in the index.
    K root_key;                         ///< The key of the root segment.
    Floating root_slope;                ///< The slope of the root segment.
    Floating root_intercept;            ///< The intercept of the root segment.
    K first_key;                        ///< The smallest element in the data.
    K last_key;                         ///< The largest element in the data.

    struct CompressedLayer {
        std::vector<K> keys;                       ///< The keys of the segment in this layer.
        PackedVector slopes_map;                   ///< The ith element is an index into slopes_table.
        std::vector<Floating> &slopes_table;       ///< A reference to the vector containing the slopes.
        int64_t intercept_offset;                  ///< An offset to make the intercepts start from 0 in the bitvector.
        CompressedBV compressed_intercepts;        ///< The compressed bitvector storing the intercepts.
        typename CompressedBV::select_1_type sel1; ///< The select1 succinct data structure on compressed_intercepts.

        template<typename IterK, typename IterI, typename IterM>
        CompressedLayer(IterK first_key, IterK last_key,
                        IterI first_intercept, IterI last_intercept,
                        IterM first_slope, IterM last_slope,
                        std::vector<Floating> &slopes_table)
            : keys(first_key, last_key), intercept_offset(*first_intercept), slopes_table(slopes_table) {
            sdsl::bit_vector intercept_bv(*(last_intercept - 1) - intercept_offset + 1);
            for (auto it = first_intercept; it != last_intercept; ++it)
                intercept_bv[*it - intercept_offset] = 1;
            compressed_intercepts = CompressedBV(intercept_bv);
            slopes_map = PackedVector(first_slope, last_slope, PackedVector::bit_length(slopes_table.size() - 1));
            sdsl::util::init_support(sel1, &compressed_intercepts);
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
            return keys.size() * sizeof(K) + slopes_map.size_in_bytes() + sdsl::size_in_bytes(compressed_intercepts);
        }
    };

    std::vector<CompressedLayer> layers;   ///< the layers composing the compressed index

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
    CompressedPGMIndex(Iterator first, Iterator last) : data_size(std::distance(first, last)) {
        assert(std::is_sorted(data.cbegin(), data.cend()));
        if (data_size > 0) {
            first_key = *first;
            last_key = *std::prev(last);
        }

        std::vector<size_t> layers_boundaries({0});
        std::vector<K> tmp_keys;
        std::vector<std::pair<Floating, Floating>> tmp_ranges;
        std::vector<std::pair<Floating, Floating>> tmp_intersections;
        tmp_keys.reserve(2048);
        tmp_ranges.reserve(2048);
        tmp_intersections.reserve(2048);
        bool needs_more_layers = true;

        // Iterative construction of the levels
        while (needs_more_layers) {
            bool first_level = layers_boundaries.size() == 1;
            size_t error = first_level ? Error : RecursiveError;
            size_t offset = first_level ? 0 : layers_boundaries[layers_boundaries.size() - 2];
            size_t start = offset;
            size_t n = first_level ? data_size : tmp_keys.size();
            auto it = first_level ? first : tmp_keys.begin() + offset;
            auto key = *it;

            // Compute segmentation for this level
            OptimalPiecewiseLinearModel<K, size_t> algorithm(error, error);
            algorithm.add_point(*it, 0);
            ++it;

            for (size_t i = offset + 1; i < n; ++i, ++it) {
                if (*it == *std::prev(it))
                    continue;

                if (!algorithm.add_point(*it, i - offset)) {
                    tmp_keys.emplace_back(key);
                    tmp_ranges.emplace_back(algorithm.get_slope_range());
                    tmp_intersections.emplace_back(algorithm.get_intersection());

                    key = *it;
                    start = i;
                    --i;
                    --it;
                }
            }

            // Last segment
            if (start < n - 1) {
                tmp_keys.emplace_back(key);
                tmp_ranges.emplace_back(algorithm.get_slope_range());
                tmp_intersections.emplace_back(algorithm.get_intersection());
            }

            needs_more_layers = tmp_keys.size() - layers_boundaries.back() > 1;
            layers_boundaries.push_back(tmp_keys.size());
        }

        // Compress the slopes
        auto[tmp_table, tmp_map, tmp_intercepts] = merge_slopes(tmp_ranges, tmp_intersections, tmp_keys);
        slopes_table = tmp_table;

        // Build layers
        auto[i_x, i_y] = tmp_intersections.back();
        root_key = tmp_keys.back();
        root_slope = 0.5 * (tmp_ranges.back().first + tmp_ranges.back().second);
        root_intercept = std::floor(i_y - (i_x - root_key) * root_slope);

        layers.reserve(layers_boundaries.size() - 2);
        for (int i = layers_boundaries.size() - 2; i > 0; --i) {
            auto l = layers_boundaries[i - 1];
            auto r = layers_boundaries[i];
            layers.emplace_back(tmp_keys.begin() + l, tmp_keys.begin() + r,
                                tmp_intercepts.begin() + l, tmp_intercepts.begin() + r,
                                tmp_map.begin() + l, tmp_map.begin() + r, slopes_table);
        }
    }

    /**
     * Returns the size in bytes of the index.
     * @return the size in bytes of the index
     */
    size_t size_in_bytes() const {
        size_t accum = 0;
        for (auto &l : layers)
            accum += l.size_in_bytes();
        return accum + slopes_table.size() * sizeof(Floating);
    }

    /**
     * Returns the approximate position of a key.
     * @param key the value to search for
     * @return a struct with the approximate position
     * @see approx_pos_t
     */
    inline ApproxPos find_approximate_position(K key) const {
        if (UNLIKELY(key < first_key))
            return {0, 0, 0};
        if (UNLIKELY(key > last_key))
            return {data_size - 1, data_size - 1, data_size - 1};

        size_t approx_pos = std::max<Floating>(0, root_intercept + root_slope * (key - root_key));
        size_t pos = 0;

        for (auto layer_i = 0; layer_i < layers.size(); ++layer_i) {
            auto &it = layers[layer_i];
            auto layer_size = it.size();
            auto lo = SUB_ERR(approx_pos, RecursiveError, layer_size);
            auto hi = ADD_ERR(approx_pos, RecursiveError + 1, layer_size);

            for (; lo <= hi && it.keys[lo] <= key; ++lo);
            pos = lo - 1;

            approx_pos = it(pos, key);
            if (pos + 1 < it.size())
                approx_pos = std::min(approx_pos, (size_t) it.get_intercept(pos + 1));
            assert(it.keys[pos] <= key);
            assert(pos + 1 >= hi || it.keys[pos + 1] > key);
        }

        auto p = layers.back()(pos, key);
        auto lo = SUB_ERR(p, Error, data_size);
        auto hi = ADD_ERR(p, Error + 1, data_size);

        return {p, lo, hi};
    }

    /**
     * Returns the number of segments in the last level of the index.
     * @return the number of segments
     */
    size_t segments_count() const {
        return layers.back().size();
    }

    /**
     * Returns the height of the index.
     * @return the height of the index
     */
    size_t height() const {
        return layers.size() + 1;
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