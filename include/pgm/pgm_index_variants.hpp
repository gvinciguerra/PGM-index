// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2020 Giorgio Vinciguerra.
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

#include <string>
#include <stdexcept>
#include "sdsl.hpp"
#include "pgm_index.hpp"

namespace pgm {

/** Computes the smallest integral value not less than x / y, where x and y must be positive integers. */
#define CEIL_UINT_DIV(x, y) ((x) / (y) + ((x) % (y) != 0))

/** Computes the number of bits needed to store x, that is, 0 if x is 0, 1 + floor(log2(x)) otherwise. */
#define BIT_WIDTH(x) ((x) == 0 ? 0 : 64 - __builtin_clzll(x))

/**
 * A simple variant of @ref BinarySearchBasedPGMIndex that builds a top-level lookup table to speed up the search on the
 * segments.
 *
 * If properly tuned and used in a system with a fast internal memory, this variant can be faster than @ref PGMIndex in
 * the average case.
 *
 * The @p TopLevelBitSize template argument allows to specify the bit-size of the memory cells in the top-level table.
 * It can be set either to a power of two or to 0. If set to 0, the bit-size of the cells will be determined dynamically
 * so that the table is bit-compressed.
 *
 * @tparam K the type of the indexed keys
 * @tparam Epsilon controls the size of the returned search range
 * @tparam TopLevelBitSize the bit-size of the cells in the top-level table, must be either 0 or a power of two
 * @tparam Floating the floating-point type to use for slopes
 */
template<typename K, size_t Epsilon = 64, uint8_t TopLevelBitSize = 32, typename Floating = float>
class InMemoryPGMIndex {
protected:
    static_assert(Epsilon > 0);
    static_assert(TopLevelBitSize == 0 || (TopLevelBitSize & (TopLevelBitSize - 1u)) == 0);

    using Segment = typename PGMIndex<K, Epsilon, 0, Floating>::Segment;

    size_t n;                                     ///< The number of elements this index was built on.
    K first_key;                                  ///< The smallest element.
    std::vector<Segment> segments;                ///< The segments composing the index.
    sdsl::int_vector<TopLevelBitSize> top_level;  ///< The structure on the segment.
    K step;

    template<typename RandomIt>
    void build_top_level(RandomIt, RandomIt last, size_t top_level_size) {
        auto log_segments = (size_t) BIT_WIDTH(segments.size());
        if constexpr (TopLevelBitSize == 0)
            top_level = sdsl::int_vector<>(top_level_size, segments.size(), log_segments);
        else {
            if (log_segments > TopLevelBitSize)
                throw std::invalid_argument("The value TopLevelBitSize=" + std::to_string(TopLevelBitSize) +
                    " is too low. Try to set it to " + std::to_string(TopLevelBitSize << 1));
            top_level = sdsl::int_vector<TopLevelBitSize>(top_level_size, segments.size(), TopLevelBitSize);
        }

        step = (size_t) CEIL_UINT_DIV(*std::prev(last), top_level_size);
        for (auto i = 0ull, k = 1ull; i < top_level_size - 1; ++i) {
            while (k < segments.size() && segments[k].key < (i + 1) * step)
                ++k;
            top_level[i] = k;
        }
    }

    /**
     * Returns the segment responsible for a given key, that is, the rightmost segment having key <= the sought key.
     * @param key the value of the element to search for
     * @return an iterator to the segment responsible for the given key
     */
    auto segment_for_key(const K &key) const {
        auto j = std::min<size_t>(key / step, top_level.size() - 1);
        auto first = segments.begin() + (key < step ? 1 : top_level[j - 1]);
        auto last = segments.begin() + top_level[j];
        return std::prev(std::upper_bound(first, last, key));
    }

public:

    static constexpr size_t epsilon_value = Epsilon;

    /**
     * Constructs an empty index.
     */
    InMemoryPGMIndex() = default;

    /**
     * Constructs the index on the given sorted vector, with the specified top level size.
     * @param data the vector of keys, must be sorted
     * @param top_level_size the number of cells allocated for the top-level table
     */
    explicit InMemoryPGMIndex(const std::vector<K> &data, size_t top_level_size)
        : InMemoryPGMIndex(data.begin(), data.end(), top_level_size) {}

    /**
     * Constructs the index on the sorted keys in the range [first, last), with the specified top level size.
     * @param first, last the range containing the sorted keys to be indexed
     * @param top_level_size the number of cells allocated for the top-level table
     */
    template<typename RandomIt>
    InMemoryPGMIndex(RandomIt first, RandomIt last, size_t top_level_size)
        : n(std::distance(first, last)),
          first_key(n ? *first : 0),
          segments(),
          top_level() {
        std::vector<size_t> sizes;
        std::vector<size_t> offsets;
        PGMIndex<K, Epsilon, 0, Floating>::build(first, last, Epsilon, 0, segments, sizes, offsets);
        build_top_level(first, last, top_level_size);
    }

    /**
     * Returns the approximate position and the range where @p key can be found.
     * @param key the value of the element to search for
     * @return a struct with the approximate position and bounds of the range
     */
    ApproxPos search(const K &key) const {
        auto k = std::max(first_key, key);
        auto it = segment_for_key(k);
        auto pos = std::min<size_t>((*it)(k), std::next(it)->intercept);
        auto lo = PGM_SUB_EPS(pos, Epsilon);
        auto hi = PGM_ADD_EPS(pos, Epsilon, n);
        return {pos, lo, hi};
    }

    /**
     * Returns the number of segments in the last level of the index.
     * @return the number of segments
     */
    size_t segments_count() const {
        return segments.size();
    }

    /**
     * Returns the number of levels of the index.
     * @return the number of levels of the index
     */
    size_t height() const {
        return 1;
    }

    /**
     * Returns the size of the index in bytes.
     * @return the size of the index in bytes
     */
    size_t size_in_bytes() const {
        return segments.size() * sizeof(Segment) + top_level.size() * top_level.width();
    }
};

/**
 * A variant of @ref BinarySearchBasedPGMIndex that builds a top-level succinct structure to speed up the search on the
 * segments.
 *
 * @tparam K the type of the indexed keys
 * @tparam Epsilon controls the size of the returned search range
 * @tparam Floating the floating-point type to use for slopes
 */
template<typename K, size_t Epsilon = 64, typename Floating = float>
class EliasFanoPGMIndex {
protected:
    static_assert(Epsilon > 0);

    using Segment = typename PGMIndex<K, Epsilon, 0, Floating>::Segment;

    struct SegmentData {
        Floating slope;    ///< The slope of the segment.
        int32_t intercept; ///< The intercept of the segment.

        SegmentData() = default;

        SegmentData(Segment &s) : slope(s.slope), intercept(s.intercept) {}

        inline size_t operator()(const K &origin, const K &k) const {
            auto pos = int64_t(slope * (k - origin)) + intercept;
            return pos > 0 ? size_t(pos) : 0ull;
        }
    };

    size_t n;                                ///< The number of elements this index was built on.
    K first_key;                             ///< The smallest segment key.
    std::vector<SegmentData> segments;       ///< The segments composing the index.
    sdsl::sd_vector<> ef;                    ///< The Elias-Fano structure on the segment.
    sdsl::sd_vector<>::rank_1_type rank;     ///< The rank1 structure.

public:

    static constexpr size_t epsilon_value = Epsilon;

    /**
     * Constructs an empty index.
     */
    EliasFanoPGMIndex() = default;

    /**
     * Constructs the index on the given sorted vector.
     * @param data the vector of keys, must be sorted
     */
    explicit EliasFanoPGMIndex(const std::vector<K> &data) : EliasFanoPGMIndex(data.begin(), data.end()) {}

    /**
     * Constructs the index on the sorted keys in the range [first, last).
     * @param first, last the range containing the sorted keys to be indexed
     */
    template<typename RandomIt>
    EliasFanoPGMIndex(RandomIt first, RandomIt last)
        : n(std::distance(first, last)),
          first_key(n ? *first : 0),
          segments(),
          ef() {
        if (n == 0)
            return;

        std::vector<Segment> tmp;
        std::vector<size_t> sizes;
        std::vector<size_t> offsets;
        PGMIndex<K, Epsilon, 0, Floating>::build(first, last, Epsilon, 0, tmp, sizes, offsets);

        segments.reserve(tmp.size());
        for (auto &x: tmp) {
            segments.push_back(x);
            x.key -= first_key;
        }

        ef = decltype(ef)(tmp.begin(), std::prev(tmp.end()));
        sdsl::util::init_support(rank, &ef);
    }

    /**
     * Returns the approximate position and the range where @p key can be found.
     * @param key the value of the element to search for
     * @return a struct with the approximate position and bounds of the range
     */
    ApproxPos search(const K &key) const {
        auto k = std::max(first_key, key);
        auto[r, origin] = rank.pred(k - first_key);
        auto pos = std::min<size_t>(segments[r](origin + first_key, k), segments[r + 1].intercept);
        auto lo = PGM_SUB_EPS(pos, Epsilon);
        auto hi = PGM_ADD_EPS(pos, Epsilon, n);
        return {pos, lo, hi};
    }

    /**
     * Returns the number of segments in the last level of the index.
     * @return the number of segments
     */
    size_t segments_count() const {
        return segments.size();
    }

    /**
     * Returns the number of levels of the index.
     * @return the number of levels of the index
     */
    size_t height() const {
        return 1;
    }

    /**
     * Returns the size of the index in bytes.
     * @return the size of the index in bytes
     */
    size_t size_in_bytes() const {
        return segments.size() * sizeof(SegmentData) + sdsl::size_in_bytes(ef);
    }
};

}