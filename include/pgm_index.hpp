// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2018 Giorgio Vinciguerra.
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

#include <limits>
#include <vector>
#include <cassert>
#include <utility>
#include <algorithm>
#include "piecewise_linear_model.hpp"

#define ADD_ERR(x, error, size) ((x) + (error) >= (size) ? (size) : (x) + (error))
#define SUB_ERR(x, error) ((x) <= (error) ? 0 : ((x) - (error)))

/**
 * A struct that stores the result of a query to a @ref PGMIndex, that is, a range [@ref lo, @ref hi)
 * centered around an approximate position @ref pos of the sought key.
 */
struct ApproxPos {
    size_t pos; ///< The approximate position of the key.
    size_t lo;  ///< The lower bound of the range where the key can be found.
    size_t hi;  ///< The upper bound of the range where the key can be found.
};

/**
 * A space-efficient index that finds the position of a key within a radius of @p Error.
 *
 * The index is constructed on a sorted sequence of keys. A query returns a struct @ref ApproxPos containing an
 * approximate position of the sought key and the bounds of the range of size 2*Error where the sought key is guaranteed
 * to be found if present. In the case of repeated keys, the index finds the position of the first occurrence of a key.
 *
 * The @p Error template parameter should be set according to the desired space-time trade-off. A smaller error value
 * makes the estimation more precise and the range smaller but at the cost of increased space usage.
 *
 * Internally the index uses a succinct piecewise linear mapping from keys to their position in the sorted order.
 * This mapping is represented as a sequence of linear models (segments) which, if @p RecursiveError is not zero, are
 * themselves recursively indexed by other piecewise linear mappings.
 *
 * @tparam K the type of the indexed elements
 * @tparam Error the maximum error allowed in the last level of the index
 * @tparam RecursiveError the maximum error allowed in the upper levels of the index
 * @tparam Floating the floating-point type to use for slopes
 */
template<typename K, size_t Error = 64, size_t RecursiveError = 16, typename Floating = double>
class PGMIndex {
    static_assert(Error > 0);
    struct Segment;

    size_t n;                           ///< The number of elements this index was built on.
    K first_key;                        ///< The smallest element.
    std::vector<Segment> segments;      ///< The segments composing the index.
    std::vector<size_t> levels_sizes;   ///< The number of segment in each level, in reverse order.
    std::vector<size_t> levels_offsets; ///< The starting position of each level in segments[], in reverse order.

    /**
     * Returns the segment responsible for a given key, that is, the rightmost segment having key <= the sought key.
     * @param key the value of the element to search for
     * @return an iterator to the segment responsible for the given key
     */
    auto segment_for_key(const K &key) const {
        if (RecursiveError == 0) {
            auto it = std::upper_bound(segments.begin(), segments.begin() + levels_sizes[0], key);
            return it == segments.begin() ? it : std::prev(it);
        }

        auto it = segments.begin() + levels_offsets.back();

        for (auto l = int(height()) - 2; l >= 0; --l) {
            auto level_begin = segments.begin() + levels_offsets[l];
            auto pos = std::min<size_t>((*it)(key), std::next(it)->intercept);
            auto lo = level_begin + SUB_ERR(pos, RecursiveError + 1);

            static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Segment);
            if constexpr (RecursiveError <= linear_search_threshold) {
                for (; std::next(lo)->key <= key; ++lo);
                it = lo;
            } else {
                auto level_size = levels_sizes[l];
                auto hi = level_begin + ADD_ERR(pos, RecursiveError + 2, level_size);
                it = std::upper_bound(lo, hi, key);
                it = it == level_begin ? it : std::prev(it);
            }
        }
        return it;
    }

public:

    static constexpr size_t error_value = Error;

    /**
     * Constructs an empty index.
     */
    PGMIndex() = default;

    /**
     * Constructs the index on the given sorted data.
     * @param data the vector of keys, must be sorted
     */
    explicit PGMIndex(const std::vector<K> &data) : PGMIndex(data.begin(), data.end()) {}

    /**
     * Constructs the index on the sorted data in the range [first, last).
     * @param first, last the range containing the sorted elements to be indexed
     */
    template<typename RandomIt>
    PGMIndex(RandomIt first, RandomIt last)
        : n(std::distance(first, last)), first_key(*first), segments(), levels_sizes(), levels_offsets() {
        assert(std::is_sorted(first, last));
        if (n == 0)
            return;

        levels_offsets.push_back(0);
        segments.reserve(n / (Error * Error));

        auto n_segments = 0ull;
        auto ignore_last = *std::prev(last) == std::numeric_limits<K>::max(); // max is reserved for padding
        auto last_n = n - ignore_last;
        last -= ignore_last;

        auto back_check = [this, last, &n_segments, &last_n]() {
            if (segments.back().slope == 0) {
                // Here, we need to ensure that keys > *(last-1) are approximated to a position == n
                segments.emplace_back(*std::prev(last) + 1, 0, last_n);
                ++n_segments;
            }
            segments.emplace_back(last_n);
        };

        // Build first level
        auto in_fun = [this, first](auto i) {
            auto x = first[i];
            if (i > 0 && i + 1u < n && x == first[i - 1] && x != first[i + 1] && x + 1 != first[i + 1])
                return std::pair<K, size_t>(x + 1, i);
            return std::pair<K, size_t>(x, i);
        };
        auto out_fun = [this](auto, auto, auto cs) { segments.emplace_back(cs); };

        n_segments = make_segmentation_par(last_n, Error, in_fun, out_fun);
        back_check();
        levels_offsets.push_back(levels_offsets.back() + n_segments + 1);
        levels_sizes.push_back(n_segments);
        last_n = n_segments;

        // Build upper levels
        while (RecursiveError && last_n > 1) {
            auto offset = levels_offsets[levels_offsets.size() - 2];
            auto in_fun_rec = [this, offset](auto i) { return std::pair<K, size_t>(segments[offset + i].key, i); };
            n_segments = make_segmentation(last_n, RecursiveError, in_fun_rec, out_fun);
            back_check();
            levels_offsets.push_back(levels_offsets.back() + n_segments + 1);
            levels_sizes.push_back(n_segments);
            last_n = n_segments;
        }

        levels_offsets.pop_back();
    }

    /**
     * Returns the approximate position of a key.
     * @param key the value of the element to search for
     * @return a struct with the approximate position
     */
    ApproxPos find_approximate_position(const K &key) const {
        auto k = std::max(first_key, key);
        auto it = segment_for_key(k);
        auto pos = std::min<size_t>((*it)(k), std::next(it)->intercept);
        auto lo = SUB_ERR(pos, Error);
        auto hi = ADD_ERR(pos, Error + 1, n);
        return {pos, lo, hi};
    }

    /**
     * Returns the number of segments in the last level of the index.
     * @return the number of segments
     */
    size_t segments_count() const {
        return segments.empty() ? 0 : levels_sizes.front();
    }

    /**
     * Returns the number of levels in the index.
     * @return the number of levels in the index
     */
    size_t height() const {
        return levels_sizes.size();
    }

    /**
     * Returns the size of the index in bytes.
     * @return the size of the index in bytes
     */
    size_t size_in_bytes() const {
        return segments.size() * sizeof(Segment);
    }
};

#pragma pack(push, 1)

template<typename K, size_t Error, size_t RecursiveError, typename Floating>
struct PGMIndex<K, Error, RecursiveError, Floating>::Segment {
    K key;             ///< The first key that the segment indexes.
    Floating slope;    ///< The slope of the segment.
    int32_t intercept; ///< The intercept of the segment.

    Segment() = default;

    Segment(size_t n) : key(std::numeric_limits<K>::max()), slope(), intercept(n) {};

    /**
     * Constructs a new segment.
     * @param key the first key that the segment indexes
     * @param slope the slope of the segment
     * @param intercept the intercept of the segment
     */
    Segment(K key, Floating slope, Floating intercept) : key(key), slope(slope), intercept(intercept) {};

    explicit Segment(const typename OptimalPiecewiseLinearModel<K, size_t>::CanonicalSegment &cs)
        : key(cs.get_first_x()) {
        auto[cs_slope, cs_intercept] = cs.get_floating_point_segment(key);
        if (cs_intercept > std::numeric_limits<decltype(intercept)>::max())
            throw std::overflow_error("Change the type of Segment::intercept to int64");
        slope = cs_slope;
        intercept = std::round(cs_intercept);
    }

    friend inline bool operator<(const Segment &s, const K &k) {
        return s.key < k;
    }

    friend inline bool operator<(const K &k, const Segment &s) {
        return k < s.key;
    }

    /**
     * Returns the approximate position of the specified key.
     * @param k the key whose position must be approximated
     * @return the approximate position of the specified key
     */
    inline size_t operator()(const K &k) const {
        auto pos = int64_t(slope * (k - key)) + intercept;
        return pos > 0 ? size_t(pos) : 0ull;
    }
};

#pragma pack(pop)

/**
 * A space-efficient index that finds the position of a sought key within a radius of @p Error. This variant uses a
 * binary search in the last level, and it should only be used when BinarySearchBasedPGMIndex::size_in_bytes() is low
 * (for example, less than the last level cache size).
 * @tparam K the type of the indexed elements
 * @tparam Error the maximum allowed error in the last level of the index
 * @tparam Floating the floating-point type to use for slopes
 */
template<typename K, size_t Error, typename Floating = double>
using BinarySearchBasedPGMIndex = PGMIndex<K, Error, 0, Floating>;
