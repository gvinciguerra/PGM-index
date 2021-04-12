// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2018 Giorgio Vinciguerra.
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

#include "pgm_index.hpp"
#include "eytzinger_array.hpp"
#include "piecewise_linear_model.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pgm {

#define PGM_SUB_EPS(x, epsilon) ((x) <= (epsilon) ? 0 : ((x) - (epsilon)))
#define PGM_ADD_EPS(x, epsilon, size) ((x) + (epsilon) + 2 >= (size) ? (size) : (x) + (epsilon) + 2)

/**
 * A space-efficient index that enables fast search operations on a sorted sequence of @c n numbers.
 *
 * A search returns a struct @ref ApproxPos containing an approximate position of the sought key in the sequence and
 * the bounds of a range of size 2*Epsilon+1 where the sought key is guaranteed to be found if present.
 * If the key is not present, the range is guaranteed to contain a key that is not less than (i.e. greater or equal to)
 * the sought key, or @c n if no such key is found.
 * In the case of repeated keys, the index finds the position of the first occurrence of a key.
 *
 * The @p Epsilon template parameter should be set according to the desired space-time trade-off. A smaller value
 * makes the estimation more precise and the range smaller but at the cost of increased space usage.
 *
 * Internally the index uses a succinct piecewise linear mapping from keys to their position in the sorted order.
 * This mapping is represented as a sequence of linear models (segments) which, if @p EpsilonRecursive is not zero, are
 * themselves recursively indexed by other piecewise linear mappings.
 *
 * @tparam K the type of the indexed keys
 * @tparam Epsilon controls the size of the returned search range
 * @tparam EpsilonRecursive controls the size of the search range in the internal structure
 * @tparam Floating the floating-point type to use for slopes
 */
template<typename K, size_t Epsilon = 64, typename Floating = float>
class PGMIndexEytzinger : public PGMIndex<K, Epsilon, 0, Floating> {
protected:
    struct IndexIterator;
    struct IdxHolder;

    using Segment = typename PGMIndex<K, Epsilon, 0,Floating>::Segment;

    EytzingerArray<IdxHolder> eytzinger_first_layer;
    /**
     * Returns the segment responsible for a given key, that is, the rightmost segment having key <= the sought key.
     * @param key the value of the element to search for
     * @return an iterator to the segment responsible for the given key
     */
    auto segment_for_key(const K &key) const {
      auto idx = eytzinger_first_layer.search(IdxHolder(key, 0));
      auto it = this->segments.begin() + ((idx == eytzinger_first_layer.size()) ? this->segments_count() : eytzinger_first_layer[idx].idx);
      return it == this->segments.begin() ? it : std::prev(it);
    }

public:
    /**
     * Constructs an empty index.
     */
    PGMIndexEytzinger() = default;

    /**
     * Constructs the index on the given sorted vector.
     * @param data the vector of keys to be indexed, must be sorted
     */
    explicit PGMIndexEytzinger(const std::vector<K> &data) : PGMIndex<K, Epsilon, 0, Floating>(data.begin(), data.end()) {}

    /**
     * Constructs the index on the sorted keys in the range [first, last).
     * @param first, last the range containing the sorted keys to be indexed
     */
    template<typename RandomIt>
    PGMIndexEytzinger(RandomIt first, RandomIt last) :
        PGMIndex<K, Epsilon, 0, Floating>(first, last), eytzinger_first_layer() {
        auto iter = IndexIterator(this->segments.begin());
        eytzinger_first_layer = EytzingerArray<IdxHolder>(iter, this->segments_count());
    }

    /**
     * Returns the size of the index in bytes.
     * @return the size of the index in bytes
     */
    size_t size_in_bytes() const { return this->segments.size() * sizeof(Segment) + this->levels_offsets.size() * sizeof(size_t) + eytzinger_first_layer.size() * sizeof(IdxHolder); }
};

#pragma pack(push, 1)

template<typename K, size_t Epsilon, typename Floating>
struct PGMIndexEytzinger<K, Epsilon, Floating>::IndexIterator {
  using T = typename std::vector<Segment>::iterator;
  T iter;
  size_t idx = 0;

 public:

  explicit IndexIterator(T iter) : iter(iter) {}

  IndexIterator & operator++() {
    ++iter;
    ++idx;
    return *this;
  }

  IdxHolder operator*() {
    return IdxHolder(iter->key, idx);
  }
};

template<typename K, size_t Epsilon, typename Floating>
struct PGMIndexEytzinger<K, Epsilon, Floating>::IdxHolder {
  K value;
  size_t idx;

  IdxHolder() = default;

  IdxHolder(K value, size_t idx) : value(value), idx(idx) {}

  IdxHolder(IdxHolder const &) = default;
  IdxHolder(IdxHolder &&) = default;
  IdxHolder & operator=(IdxHolder const &) noexcept = default;
  IdxHolder & operator=(IdxHolder &&) noexcept = default;

  bool operator<(IdxHolder const &val) {
    return value < val.value;
  }
  bool operator<=(IdxHolder const &val) {
    return value <= val.value;
  }
};

#pragma pack(pop)

}