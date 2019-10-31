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

#include <list>
#include <tuple>
#include <cmath>
#include <vector>
#include <cstring>
#include <cassert>
#include <algorithm>
#include "piecewise_linear_model.hpp"

#define ADD_ERR(x, error, size) ((x) + (error) >= (size) ? (size) - 1 : (x) + (error))
#define SUB_ERR(x, error, size) ((x) <= (error) ? 0 : ((x) - (error)))
#define BIN_SEARCH_THRESHOLD 256

#define LIKELY(x)   __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

/**
 * A struct that stores the result of a query to a PGMIndex.
 */
struct ApproxPos {
    size_t pos; ///< The approximate position.
    size_t lo;  ///< The lower bound of the range of size no more than 2*error where key can be found.
    size_t hi;  ///< The upper bound of the range of size no more than 2*error where key can be found.
};

/**
 * A struct that stores the parameters (slope and intercept) of a segment.
 * @tparam Floating the floating-point type of the segment's parameters
 */
template<typename Floating>
struct SegmentData {
    static_assert(std::is_floating_point<Floating>());
    Floating slope;     ///< The slope of the segment.
    Floating intercept; ///< The intercept of the segment.

    SegmentData(Floating slope, Floating intercept) : slope(slope), intercept(intercept) {};

    template<typename K>
    inline size_t operator()(K k) const {
        Floating pos = slope * k + intercept;
        return pos > Floating(0) ? pos : 0ul;
    }
};

/**
 * A struct that stores a segment.
 * @tparam K the type of the elements that the segment indexes
 * @tparam Floating the floating-point type of the segment's parameters
 */
template<typename K, typename Floating>
struct Segment {
    static_assert(std::is_floating_point<Floating>());
    K key;              ///< The first key that the segment indexes.
    Floating slope;     ///< The slope of the segment.
    Floating intercept; ///< The intercept of the segment.

    Segment() = default;

    /**
     * Constructs a new segment.
     * @param key the first key that the segment indexes
     * @param slope the slope of the segment
     * @param intercept the intercept of the segment
     */
    Segment(K key, Floating slope, Floating intercept) : key(key), slope(slope), intercept(intercept) {};

    friend inline bool operator<(const Segment &s, const K k) {
        return s.key < k;
    }

    friend inline bool operator<(const Segment &s1, const Segment &s2) {
        return s1.key < s2.key;
    }

    /**
     * Returns the approximate position of the specified key.
     * @param k the key whose position must be approximated
     * @return the approximate position of the specified key
     */
    inline size_t operator()(K k) const {
        Floating pos = slope * (k - key) + intercept;
        return pos > Floating(0) ? pos : 0ul;
    }
};

/**
 * Stores the result of a segmentation, that is, a piecewise linear model mapping keys to locations.
 * @tparam K the type of the elements upon which the segmentation is computed
 * @tparam Error the maximum allowed error of the segmentation
 * @tparam Floating the floating-point type used for slopes and intercepts
 */
template<typename K, size_t Error, typename Floating = double>
struct Segmentation {
    using floating_type = Floating;
    using segment_type = Segment<K, Floating>;
    using segment_data_type = SegmentData<Floating>;
    using key_type = K;

    static const size_t error = Error;  ///< The maximum error of the segmentation
    size_t data_size = 0;               ///< The size of the data upon which the segmentation was computed.
    std::vector<segment_type> segments; ///< The vector of segments.

    Segmentation() = default;

    /**
    * Builds a piecewise linear model with the given data.
    * @param data the vector of keys, must be sorted
    * @return the vector of segments approximating the data
    */
    explicit Segmentation(const std::vector<K> &data)
        : data_size(data.size()), segments(build_segments(data, error)) {}

    /**
     * Computes a piecewise linear model with the given data.
     * @param data the vector of keys, must be sorted
     * @param error the maximum allowed error of the piecewise linear model
     * @return the piecewise linear model as a vector of segments
     */
    static std::vector<segment_type> build_segments(const std::vector<K> &data, size_t error) {
        if (data.empty())
            return {};

        std::vector<segment_type> segments;
        segments.reserve(8192);

        OptimalPiecewiseLinearModel<K, size_t> algorithm(error, error);
        algorithm.add_point(data[0], 0);
        size_t start = 0;

        for (size_t i = 1; i < data.size(); ++i) {
            if (data[i] == data[i - 1])
                continue;

            if (!algorithm.add_point(data[i], i)) {
                auto key = data[start];
                auto intercept = algorithm.get_intercept(key);
                auto[min_slope, max_slope] = algorithm.get_slope_range();
                auto slope = 0.5 * (min_slope + max_slope);
                segments.emplace_back(key, slope, intercept);
                start = i;
                --i;
            }
        }

        // Last segment
        if (start < data.size() - 1) {
            auto key = data[start];
            auto intercept = algorithm.get_intercept(key);
            auto[min_slope, max_slope] = algorithm.get_slope_range();
            auto slope = 0.5 * (min_slope + max_slope);
            segments.emplace_back(key, slope, intercept);
        }

        return segments;
    }

};

/**
 * An indexing strategy for PGMIndex that stores the segment in an index built by repeating the segmentation process
 * until there is only one segment left.
 * @tparam RecursiveError the error of the repeated segmentation
 */
template<typename SegmentationType, size_t RecursiveError>
class RecursiveStrategy {
    using K = typename SegmentationType::key_type;
    using Floating = typename SegmentationType::floating_type;
    using segment_type = typename SegmentationType::segment_type;
    using segment_data_type = SegmentData<Floating>;
    segment_type root;
    size_t root_limit;

    struct Layer {
        std::vector<K> segments_keys;
        std::vector<segment_data_type> segments_data;

        inline size_t size() const {
            return segments_keys.size();
        }

        template<typename S>
        explicit Layer(const S &segmentation) {
            segments_keys.reserve(segmentation.segments.size());
            segments_data.reserve(segmentation.segments.size() + 1);
            for (auto &s : segmentation.segments) {
                segments_keys.push_back(s.key);
                segments_data.emplace_back(s.slope, s.intercept);
            }
            segments_data.emplace_back(0, segmentation.data_size);
        }
    };

    std::vector<Layer> layers;

protected:

    RecursiveStrategy() = default;

    RecursiveStrategy(const std::vector<K> &data, const SegmentationType &segmentation) : layers() {
        if (segmentation.segments.empty()) {
            root = segment_type(0, 0, 0);
            root_limit = 1;
            return;
        }

        if (segmentation.segments.size() == 1) {
            root = segment_type(segmentation.segments[0]);
            root_limit = data.size();
            return;
        }

        std::list<Layer> tmp;
        tmp.emplace_front(segmentation);
        Segmentation<K, RecursiveError, Floating> l;

        while (tmp.front().size() > 1) {
            l = Segmentation<K, RecursiveError, Floating>(tmp.front().segments_keys);
            tmp.emplace_front(l);
        }

        root = l.segments[0];
        root_limit = l.data_size;
        layers = {std::make_move_iterator(std::next(tmp.begin())), std::make_move_iterator(tmp.end())};
    }

public:

    /**
     * Finds the last-level segment responsible for the given key by searching in the recursive PGM-index structure.
     * @param key the value to search for
     * @return the segment responsible for the given key
     */
    inline const segment_type find_segment_for_key(K key) const {
        auto slope = root.slope;
        auto intercept = root.intercept;
        auto node_key = root.key;
        size_t approx_pos = std::min(root(key), root_limit);
        size_t pos = 0;

        for (auto &it : layers) {
            auto layer_size = it.size();

            auto lo = SUB_ERR(approx_pos, RecursiveError, layer_size);
            auto hi = ADD_ERR(approx_pos, RecursiveError + 1, layer_size);
            assert(it.segments_keys[lo] <= key);
            assert(key <= it.segments_keys[hi] || key > it.segments_keys[layer_size - 1]);

            if constexpr (RecursiveError >= BIN_SEARCH_THRESHOLD) { // use binary search for large "pages"
                auto layer_begin = it.segments_keys.cbegin();
                auto lo_it = layer_begin + lo;
                auto hi_it = layer_begin + hi;
                auto pos_it = std::lower_bound(lo_it, hi_it + 1, key);
                pos = (size_t) std::distance(layer_begin, pos_it);
                if (layer_begin[pos] > key && pos != lo)
                    pos--;
            } else {
                for (; lo <= hi && it.segments_keys[lo] <= key; ++lo);
                pos = lo - 1;
            }

            node_key = it.segments_keys[pos];
            approx_pos = it.segments_data[pos](key - node_key);
            slope = it.segments_data[pos].slope;
            intercept = it.segments_data[pos].intercept;
            if (pos + 1 <= it.size())
                approx_pos = std::min(approx_pos, (size_t) it.segments_data[pos + 1].intercept);

            assert(node_key <= key);
            assert(pos + 1 >= hi || it.segments_keys[pos + 1] > key);
        }

        return {node_key, slope, intercept};
    }

    /**
     * Returns the size in bytes of the recursive indexing data structure, including the last-level segments.
     * @return the size in bytes of the indexing data structure
     */
    size_t size_in_bytes() const {
        auto total = 1;
        for (auto &l : layers)
            total += l.size();
        return total * sizeof(segment_type);
    }

    /**
     * Returns the number of segments in the last layer of the approximate index.
     * @return the number of segments
     */
    size_t segments_count() const {
        return layers.empty() ? 1 : layers.back().size();
    }

    /**
     * Returns the number of layers of the recursive segmentation.
     * @return the number of layers of the recursive segmentation
     */
    size_t height() const {
        return 1 + layers.size();
    }

};

/**
 * An indexing strategy for PGMIndex that stores the segments in a sorted vector and, at query time, performs a binary
 * search on it.
 */
template<typename SegmentationType>
class BinarySearchStrategy {
    using K = typename SegmentationType::key_type;
    using segment_type = typename SegmentationType::segment_type;

    std::vector<segment_type> segments;

protected:

    BinarySearchStrategy() = default;

    BinarySearchStrategy(const std::vector<K> &data, const SegmentationType &segmentation)
        : segments(segmentation.segments) {}

public:

    /**
      * Finds the last-level segment responsible for the given key via a binary search on the array of segments.
      * @param key the value to search for
      * @return the segment responsible for the given key
      */
    inline const segment_type find_segment_for_key(K key) const {
        auto pos = std::lower_bound(segments.cbegin(), segments.cend(), key);
        auto i = std::distance(segments.cbegin(), pos) + (pos == segments.cbegin() || key == pos->key ? 0 : -1);
        return segments[i];
    }

    /**
     * Returns the size of the vector storing the segments in bytes.
     * @return the size of the vector storing the segments in bytes
     */
    size_t size_in_bytes() const {
        return segments.size() * sizeof(segment_type);
    }

    /**
     * Returns the number of segments in the approximate index.
     * @return the number of segments
     */
    size_t segments_count() const {
        return segments.size();
    }

};

/**
 * An indexing strategy for PGMIndex that stores the segments in an implicit tree.
 * @tparam NodeSize the size of a node of the tree in bytes
 */
template<typename SegmentationType, size_t NodeSize>
class TreeStrategy {
    using K = typename SegmentationType::key_type;
    using segment_type = typename SegmentationType::segment_type;
    using segment_data_type = typename SegmentationType::segment_data_type;

    size_t tree_height;
    size_t half_marker;
    size_t n_internal_nodes;
    std::vector<K> tree;
    std::vector<K> leaves;
    std::vector<segment_data_type> segments_data;
    const size_t slots_per_node = NodeSize / sizeof(K);

    template<class Iterator>
    inline const segment_type find_in_leaves(Iterator lo, Iterator hi, K key) const {
        for (; lo != hi && *lo < key; ++lo);
        long offset = std::distance(leaves.cbegin(), lo);
        auto idx = offset == 0 || key == *lo ? offset : offset - 1;
        return {leaves[idx], segments_data[idx].slope, segments_data[idx].intercept};
    }

protected:

    TreeStrategy() = default;

    TreeStrategy(const std::vector<K> &data, const SegmentationType &segmentation) : leaves() {
        const std::vector<segment_type> &segments = segmentation.segments;

        // store segments and keys of segments separately: the latter are accessed more frequently and it is
        // better to pack them together
        leaves.reserve(segments.size());
        segments_data.reserve(segments.size());
        for (auto &s : segments) {
            leaves.push_back(s.key);
            segments_data.emplace_back(s.slope, s.intercept);
        }

        const auto n = leaves.size();
        const auto leaf_nodes = ceil(n / float(slots_per_node));
        tree_height = size_t(ceil(log(leaf_nodes) / log(slots_per_node + 1)));
        const auto expp = size_t(pow(slots_per_node + 1, tree_height));
        const auto last_internal_node = size_t((expp - leaf_nodes) / slots_per_node);
        n_internal_nodes = (size_t) ((expp - 1) / slots_per_node) - last_internal_node;
        tree = std::vector<K>(n_internal_nodes * slots_per_node);
        half_marker = (expp - 1) / slots_per_node;

        auto i = tree.size();
        while (i-- > 0) {
            auto node = i / slots_per_node;
            auto child = node * (slots_per_node + 1) + 1 + i % slots_per_node;
            while (child < n_internal_nodes) // follow rightmost branch
                child = child * (slots_per_node + 1) + slots_per_node + 1;

            // child is a leaf -> map it to an index in the tree
            long diff = (child - half_marker) * slots_per_node;
            if (diff < 0)
                tree[i] = leaves[diff + n + slots_per_node - 1];
            else if (diff + slots_per_node - 1 < n - last_internal_node * slots_per_node)
                tree[i] = leaves[diff + slots_per_node - 1];
            else
                // special case: fill ancestor of the last leaf node with the
                // last element of (the biggest in) the first half of the tree
                tree[i] = leaves[n - last_internal_node * slots_per_node - 1];
        }
    }

public:

    /**
      * Finds the last-level segment responsible for the given key by traversing the tree indexing the keys of the
      * segments.
      * @param key the value to search for
      * @return the segment responsible for the given key
      */
    inline const segment_type find_segment_for_key(K key) const {
        if (n_internal_nodes == 0)
            return find_in_leaves(leaves.cbegin(), leaves.cend(), key);

        size_t child = 0;
        while (child < n_internal_nodes) {
            auto index_in_tree = child * slots_per_node;
            if (NodeSize >= BIN_SEARCH_THRESHOLD) { // use binary search for large pages and scan for smaller ones
                auto lo = tree.cbegin() + index_in_tree;
                auto hi = std::min(tree.cend(), lo + slots_per_node + 1);
                auto pos = std::lower_bound(lo, hi, key);
                if (pos == hi)
                    --pos;
                child = child * (slots_per_node + 1) + 1 + std::distance(lo, pos);
            } else {
                size_t lo;
                for (lo = 0; lo < slots_per_node && tree[index_in_tree + lo] < key; ++lo);
                child = child * (slots_per_node + 1) + 1 + lo;
            }
        }

        long diff = (long(child) - long(half_marker)) * slots_per_node;
        if (diff < 0)
            diff += leaves.size();
        assert(diff >= 0);

        auto lo = leaves.cbegin() + std::min(leaves.size(), size_t(diff));
        auto hi = leaves.cbegin() + std::min(leaves.size(), diff + slots_per_node);
        return find_in_leaves(lo, hi, key);
    }

    /**
     * Returns the overall size in bytes of the segments and the implicit tree that indexes them.
     * @return the overall size in bytes of the segments and the tree
     */
    size_t size_in_bytes() const {
        return (tree.size() + leaves.size()) * sizeof(K) + segments_data.size() * sizeof(segment_data_type);
    }

    /**
     * Returns the number of segments in the approximate index.
     * @return the number of segments
     */
    size_t segments_count() const {
        return segments_data.size();
    }

    /**
     * Returns the height of the tree.
     * @return the height of the tree
     */
    size_t height() const {
        return tree_height;
    }

};

/**
 * A space-efficient index that finds the position of a sought key within a radius of @p Error.
 *
 * Internally it uses a piecewise linear mapping from keys to their position in the sorted order. This mapping is
 * represented as a sequence of segment_t, which is itself indexed in one of three ways (collectively referred to as
 * @p IndexingStrategy): a binary search over the sequence of segments, an implicit multiway tree, a recursive
 * construction of the piecewise linear mapping (the default one). The three choices correspond to BinarySearchStrategy,
 * TreeStrategy and RecursiveStrategy, respectively.
 *
 * The type aliases @c BinarySearchBasedPGMIndex, @c TreeBasedPGMIndex and @c RecursivePGMIndex are provided for
 * convenience. For example, the following two lines are equivalent:
 *
 * ```cpp
 * PGMIndex<int, 64, RecursiveStrategy<segmentation_t<int, Error>, 16>> index(data);
 * RecursivePGMIndex<int, 64, 16> index(data);
 * ```
 *
 * @tparam K the type of the indexed elements
 * @tparam Error the maximum allowed error in the last level of the index
 * @tparam IndexingStrategy the strategy that searches the segments indexing the data
 * @tparam Floating the floating-point type to use for slopes and intercept
 */
template<typename K, size_t Error,
    typename IndexingStrategy = RecursiveStrategy<Segmentation<K, Error>, 16>,
    typename Floating = double>
class PGMIndex : public IndexingStrategy {
    using segmentation_type = Segmentation<K, Error, Floating>;

    size_t data_size; ///< The number of elements in the data.
    K first;          ///< The smallest element in the data.
    K last;           ///< The largest element in the data.

public:

    /**
     * Builds an empty index.
     */
    PGMIndex() = default;

    /**
     * Builds the index on the given data.
     * @param data the vector of keys, must be sorted
     */
    explicit PGMIndex(const std::vector<K> &data)
        : IndexingStrategy(data, segmentation_type(data)), data_size(data.size()), first(), last() {
        if (!data.empty()) {
            first = data.front();
            last = data.back();
        }
        assert(std::is_sorted(data.cbegin(), data.cend()));
    }

    /**
     * Builds the index with a precomputed segmentation.
     * @param data the vector of keys, must be sorted
     * @param segmentation the precomputed segmentation
     */
    explicit PGMIndex(const std::vector<K> &data, const segmentation_type &segmentation)
        : IndexingStrategy(data, segmentation), data_size(data.size()), first(), last() {
        if (!data.empty()) {
            first = data.front();
            last = data.back();
        }
        assert(std::is_sorted(data.cbegin(), data.cend()));
        assert(segmentation.data_size == data_size);
    }

    /**
     * Returns the approximate position of a key.
     * @param key the value to search for
     * @return a struct with the approximate position
     * @see approx_pos_t
     */
    inline ApproxPos find_approximate_position(K key) const {
        if (UNLIKELY(key < first))
            return {0, 0, 0};
        if (UNLIKELY(key > last))
            return {data_size - 1, data_size - 1, data_size - 1};

        auto segment = this->find_segment_for_key(key);
        auto pos = segment(key);
        auto lo = SUB_ERR(pos, Error, data_size);
        auto hi = ADD_ERR(pos, Error, data_size);
        if (UNLIKELY(pos > hi))
            pos = hi;

        return {pos, lo, hi};
    }

};

/**
 * A space-efficient index that finds the position of a sought key within a radius of @p Error. This variant uses a
 * recursive structure to route query keys to the last level of the index.
 * @tparam K the type of the indexed elements
 * @tparam Error the maximum allowed error in the last level of the index
 * @tparam RecursiveError the maximum allowed error in the upper levels of the index
 * @tparam Floating the floating-point type used for slopes and intercepts
 */
template<typename K, size_t Error, size_t RecursiveError, typename Floating = double>
using RecursivePGMIndex = PGMIndex<K, Error, RecursiveStrategy<Segmentation<K, Error, Floating>, RecursiveError>,
                                   Floating>;

/**
 * A space-efficient index that finds the position of a sought key within a radius of @p Error. This variant uses an
 * implicit multiway tree to route query keys to the last level of the index.
 * @tparam K the type of the indexed elements
 * @tparam Error the maximum allowed error in the last level of the index
 * @tparam NodeSize the size of a node of the tree in bytes
 * @tparam Floating the floating-point type used for slopes and intercepts
 */
template<typename K, size_t Error, size_t NodeSize, typename Floating = double>
using TreeBasedPGMIndex = PGMIndex<K, Error, TreeStrategy<Segmentation<K, Error, Floating>, NodeSize>, Floating>;

/**
 * A space-efficient index that finds the position of a sought key within a radius of @p Error. This variant uses a
 * binary search in the last level, and it should only be used when BinarySearchBasedPGMIndex::size_in_bytes() is low
 * (for example, less than the last level cache size).
 * @tparam K the type of the indexed elements
 * @tparam Error the maximum allowed error in the last level of the index
 * @tparam Floating the floating-point type used for slopes and intercepts
 */
template<typename K, size_t Error, typename Floating = double>
using BinarySearchBasedPGMIndex = PGMIndex<K, Error, BinarySearchStrategy<Segmentation<K, Error, Floating>>, Floating>;
