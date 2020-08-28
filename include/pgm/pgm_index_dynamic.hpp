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

#include <queue>
#include <limits>
#include <vector>
#include <utility>
#include <numeric>
#include <cassert>
#include <iterator>
#include <algorithm>
#include <unordered_set>
#include "pgm_index.hpp"

namespace pgm {

/**
 * A sorted associative container that contains key-value pairs with unique keys.
 * @tparam K the type of a key
 * @tparam V the type of a value
 * @tparam PGMType the type of @ref PGMIndex to use in the container
 * @tparam MinIndexedLevel the minimum level (of size 2^MinIndexedLevel) on which a @p PGMType index is constructed
 */
template<typename K, typename V, typename PGMType = PGMIndex<K>, uint8_t MinIndexedLevel = 18>
class DynamicPGMIndex {
    class ItemA;
    class ItemB;
    class Iterator;

    template<typename T, typename A=std::allocator<T>>
    class DefaultInitAllocator;

    constexpr static uint8_t min_level = 6; ///< 2^min_level-1 is the size of the sorted buffer for new items.
    constexpr static uint8_t fully_allocated_levels = std::max(15, min_level + 1);
    constexpr static size_t buffer_max_size = (1ull << (min_level + 1)) - 1;

    static_assert(min_level < MinIndexedLevel);
    static_assert(fully_allocated_levels > min_level);
    static_assert(2 * PGMType::epsilon_value < 1ul << MinIndexedLevel);

    using Item = std::conditional_t<std::is_pointer_v<V>, ItemA, ItemB>;
    using Level = std::vector<Item, DefaultInitAllocator<Item>>;

    uint8_t used_levels;       ///< Equal to 1 + last level whose size is greater than 0, or = min_level if no data.
    std::vector<Level> levels; ///< (i-min_level)th element is the data array on the ith level.
    std::vector<PGMType> pgms; ///< (i-MinIndexedLevel)th element is the index on the ith level.

    const Level &level(uint8_t level) const { return levels[level - min_level]; }
    const PGMType &pgm(uint8_t level) const { return pgms[level - MinIndexedLevel]; }
    Level &level(uint8_t level) { return levels[level - min_level]; }
    PGMType &pgm(uint8_t level) { return pgms[level - MinIndexedLevel]; }
    static uint8_t ceil_log2(size_t n) { return n <= 1 ? 0 : sizeof(unsigned long long) * 8 - __builtin_clzll(n - 1); }

    template<bool SkipDeleted, typename In1, typename In2, typename OutIterator>
    static OutIterator merge(In1 first1, In1 last1, In2 first2, In2 last2, OutIterator result) {
        while (first1 != last1 && first2 != last2) {
            if (*first2 < *first1) {
                *result = *first2;
                ++first2;
                ++result;
            } else if (*first1 < *first2) {
                *result = *first1;
                ++first1;
                ++result;
            } else if (SkipDeleted && first1->deleted()) {
                ++first1;
                ++first2;
            } else {
                *result = *first1;
                ++first1;
                ++first2;
                ++result;
            }
        }
        return std::copy(first2, last2, std::copy(first1, last1, result));
    }

    void pairwise_merge(const Item &new_item,
                        uint8_t up_to_level,
                        size_t size_hint,
                        typename Level::iterator insertion_point) {
        auto target = up_to_level + 1;
        auto actual_size = size_t(1) << (1 + min_level);
        assert((1ull << target) - level(target).size() >= actual_size);

        Level tmp_a(size_hint);
        Level tmp_b(size_hint + level(target).size());

        // Insert new_item in sorted order in the first level
        const auto tmp1 = tmp_a.begin();
        const auto tmp2 = tmp_b.begin();
        const auto mod = (up_to_level - min_level) % 2;

        auto it = std::move(levels[0].begin(), insertion_point, mod ? tmp1 : tmp2);
        *it = new_item;
        ++it;
        std::move(insertion_point, levels[0].end(), it);

        // Merge subsequent levels
        uint8_t limit = level(target).empty() ? up_to_level : up_to_level + 1;
        for (uint8_t i = 1 + min_level; i <= limit; ++i) {
            bool alternate = (i - min_level) % 2 == mod;
            auto tmp_it = alternate ? tmp1 : tmp2;
            auto out_it = alternate ? tmp2 : tmp1;
            decltype(out_it) out_end;

            bool can_delete_permanently = i == used_levels - 1;
            if (can_delete_permanently)
                out_end = merge<true>(tmp_it, tmp_it + actual_size, level(i).begin(), level(i).end(), out_it);
            else
                out_end = merge<false>(tmp_it, tmp_it + actual_size, level(i).begin(), level(i).end(), out_it);
            actual_size = std::distance(out_it, out_end);

            // Empty this level and the corresponding index
            level(i).clear();
            if (i >= fully_allocated_levels)
                level(i).shrink_to_fit();
            if (i >= MinIndexedLevel)
                pgm(i) = PGMType();
        }

        tmp_b.resize(actual_size);
        target = std::max<uint8_t>(ceil_log2(actual_size), min_level + 1);
        levels[0].resize(0);
        levels[target - min_level] = std::move(tmp_b);

        // Rebuild index, if needed
        if (target >= MinIndexedLevel)
            pgm(target) = PGMType(level(target).begin(), level(target).end());
    }

    void insert(const Item &new_item) {
        auto insertion_point = std::lower_bound(levels[0].begin(), levels[0].end(), new_item);
        if (insertion_point != levels[0].end() && *insertion_point == new_item) {
            *insertion_point = new_item;
            return;
        }

        if (levels[0].size() < buffer_max_size) {
            levels[0].insert(insertion_point, new_item);
            used_levels = used_levels == min_level ? min_level + 1 : used_levels;
            return;
        }

        size_t slots_required = buffer_max_size + 1;
        uint8_t i;
        for (i = min_level + 1; i < used_levels; ++i) {
            size_t slots_left_in_level = (1ull << i) - level(i).size();
            if (slots_required <= slots_left_in_level)
                break;
            slots_required += level(i).size();
        }

        auto insertion_idx = std::distance(levels[0].begin(), insertion_point);
        auto need_new_level = i == used_levels;
        if (need_new_level) {
            ++used_levels;
            levels.emplace_back();
            if (i - MinIndexedLevel >= int(pgms.size()))
                pgms.emplace_back();
        }

        pairwise_merge(new_item, i - 1, slots_required, levels[0].begin() + insertion_idx);
    }

public:

    using key_type = K;
    using mapped_type = V;
    using value_type = Item;
    using size_type = size_t;
    using iterator = Iterator;

    /**
     * Constructs an empty container.
     */
    DynamicPGMIndex() : used_levels(min_level), levels(32 - min_level), pgms() {
        level(min_level).reserve(buffer_max_size);
        for (uint8_t i = min_level + 1; i < fully_allocated_levels; ++i)
            level(i).reserve(1ull << i);
    }

    /**
     * Constructs the container on the sorted data in the range [first, last).
     * @tparam Iterator
     * @param first, last the range containing the sorted elements to be indexed
     */
    template<typename Iterator>
    DynamicPGMIndex(Iterator first, Iterator last) {
        size_t n = std::distance(first, last);
        used_levels = std::max<uint8_t>(ceil_log2(n), min_level) + 1;
        levels = decltype(levels)(std::max<uint8_t>(used_levels, 32) - min_level + 1);

        level(min_level).reserve(buffer_max_size);
        for (uint8_t i = min_level + 1; i < fully_allocated_levels; ++i)
            level(i).reserve(1ull << i);

        if (n == 0) {
            used_levels = min_level;
            return;
        }

        // Copy only the first of each group of pairs with same key value
        auto &target = level(used_levels - 1);
        target.resize(n);
        auto out = target.begin();
        *out++ = Item(first->first, first->second);
        while (++first != last) {
            if (first->first < std::prev(out)->first)
                throw std::invalid_argument("Range is not sorted");
            if (first->first != std::prev(out)->first)
                *out++ = Item(first->first, first->second);
        }
        target.resize(std::distance(target.begin(), out));

        if (used_levels - 1 >= MinIndexedLevel) {
            pgms = decltype(pgms)(used_levels - MinIndexedLevel);
            pgm(used_levels - 1) = PGMType(target.begin(), target.end());
        }
    }

    /**
     * Inserts an element into the container if @p key does not exists in the container. If @p key already exists, the
     * corresponding value is updated with @p value.
     * @param key element key to insert or update
     * @param value element value to insert
     */
    void insert_or_assign(const K &key, const V &value) { insert(Item(key, value)); }

    /**
     * Removes the specified element from the container.
     * @param key key value of the element to remove
     */
    void erase(const K &key) { insert(Item(key)); }

    /**
     * Finds an element with key equivalent to @p key.
     * @param key key value of the element to search for
     * @return an iterator to an element with key equivalent to @p key. If no such element is found, end() is returned
     */
    iterator find(const K &key) const {
        for (auto i = min_level; i < used_levels; ++i) {
            if (level(i).empty())
                continue;

            auto first = level(i).begin();
            auto last = level(i).end();
            if (i >= MinIndexedLevel) {
                auto range = pgm(i).search(key);
                first = level(i).begin() + range.lo;
                last = level(i).begin() + range.hi;
            }

            auto it = std::lower_bound(first, last, key);
            if (it != level(i).end() && it->first == key)
                return it->deleted() ? end() : iterator(this, it);
        }

        return end();
    }

    /**
     * Returns an iterator pointing to the first element that is not less than (i.e. greater or equal to) @p key.
     * @param key key value to compare the elements to
     * @return an iterator to an element with key not less than @p key. If no such element is found, end() is returned
     */
    iterator lower_bound(const K &key) const {
        typename Level::const_iterator lb;
        auto lb_set = false;
        std::unordered_set<K> deleted;

        for (auto i = min_level; i < used_levels; ++i) {
            if (level(i).empty())
                continue;

            auto first = level(i).begin();
            auto last = level(i).end();
            if (i >= MinIndexedLevel) {
                auto range = pgm(i).search(key);
                first = level(i).begin() + range.lo;
                last = level(i).begin() + range.hi;
            }

            auto it = std::lower_bound(first, last, key);
            for (; it != level(i).end() && it->deleted(); ++it)
                deleted.emplace(it->first);

            if (it != level(i).end() && it->first >= key && (!lb_set || it->first < lb->first)
                && deleted.find(it->first) == deleted.end()) {
                lb = it;
                lb_set = true;
            }
        }

        if (lb_set)
            return iterator(this, lb);
        return end();
    }

    /**
     * Checks if the container has no elements, i.e. whether begin() == end().
     * @return true if the container is empty, false otherwise
     */
    bool empty() const { return begin() == end(); }

    /**
     * Returns an iterator to the beginning.
     * @return an iterator to the beginning
     */
    iterator begin() const { return lower_bound(std::numeric_limits<K>::min()); }

    /**
     * Returns an iterator to the end.
     * @return an iterator to the end
     */
    iterator end() const { return iterator(this, levels.back().end()); }

    /**
     * Returns the number of elements with key that compares equal to the specified argument key, which is either 1
     * or 0 since this container does not allow duplicates.
     * @param key key value of the elements to count
     * @return number of elements with the given key, which is either 1 or 0.
     */
    size_t count(const K &key) const { return find(key) == end() ? 0 : 1; }

    /**
     * Returns the number of elements in the container.
     * @return the number of elements in the container
     */
    size_t size() const {
        // TODO: scanning the levels and using a hash table for the encountered keys could be more time efficient
        return std::distance(begin(), end());
    }

    /**
     * Returns the size of the container in bytes.
     * @return the size of the container in bytes
     */
    size_t size_in_bytes() const {
        size_t bytes = levels.size() * sizeof(Level);
        for (auto &l: levels)
            bytes += l.size() * sizeof(Item);
        return index_size_in_bytes() + bytes;
    }

    /**
     * Returns the size of the index used in this container in bytes.
     * @return the size of the index used in this container in bytes
     */
    size_t index_size_in_bytes() const {
        size_t bytes = 0;
        for (auto &p: pgms)
            bytes += p.size_in_bytes();
        return bytes;
    }

private:

    // from: https://stackoverflow.com/a/21028912
    template<typename T, typename A>
    class DefaultInitAllocator : public std::allocator<T> {
        typedef std::allocator_traits<A> a_t;

    public:
        template<typename U>
        struct rebind {
            using other = DefaultInitAllocator<U, typename a_t::template rebind_alloc<U> >;
        };

        using A::A;

        template<typename U>
        void construct(U *ptr) noexcept(std::is_nothrow_default_constructible_v<U>) {
            ::new(static_cast<void *>(ptr)) U;
        }

        template<typename U, typename...Args>
        void construct(U *ptr, Args &&... args) {
            a_t::construct(static_cast<A &>(*this), ptr, std::forward<Args>(args)...);
        }
    };
};

template<typename K, typename V, typename PGMType, uint8_t MinIndexedLevel>
class DynamicPGMIndex<K, V, PGMType, MinIndexedLevel>::Iterator {
    friend class DynamicPGMIndex;

    using internal_iterator = typename Level::const_iterator;
    using queue_pair = std::pair<internal_iterator, int16_t>;
    using dynamic_pgm_type = DynamicPGMIndex<K, V, PGMType, MinIndexedLevel>;

    static bool queue_cmp(const queue_pair &e1, const queue_pair &e2) {
        return *e1.first > *e2.first || (*e1.first == *e2.first && e1.second > e2.second);
    }

    const dynamic_pgm_type *super;
    internal_iterator it;
    bool initialized;
    std::priority_queue<queue_pair, std::vector<queue_pair>, decltype(&queue_cmp)> queue;

    void lazy_initialize_queue() {
        if (initialized)
            return;

        std::vector<queue_pair> initial_pairs{};
        initial_pairs.reserve(super->used_levels - super->min_level);

        // For each level create and position an iterator to the first key > it->first
        for (uint8_t i = super->min_level; i < super->used_levels; ++i) {
            auto &level = super->level(i);
            if (level.empty())
                continue;

            size_t lo = 0;
            size_t hi = level.size();
            if (i >= MinIndexedLevel) {
                auto range = super->pgm(i).search(it->first);
                lo = range.lo;
                hi = range.hi;
            }

            auto pos = std::upper_bound(level.begin() + lo, level.begin() + hi, *it);
            if (pos != level.end())
                initial_pairs.emplace_back(pos, i);
        }

        queue = decltype(queue)(&queue_cmp, initial_pairs);
        initialized = true;
    }

    void advance() {
        if (queue.empty()) {
            *this = super->end();
            return;
        }

        auto queue_step = [&] {
            auto[level_it, level_idx] = queue.top();
            queue.pop();
            if (std::next(level_it) != super->level(level_idx).end())
                queue.emplace(std::next(level_it), level_idx);
            return level_it;
        };

        internal_iterator tmp_it;
        do {
            tmp_it = queue_step();
            while (!queue.empty() && queue.top().first->first == tmp_it->first)
                queue_step();
        } while (!queue.empty() && tmp_it->deleted());

        if (tmp_it->deleted())
            *this = super->end();
        else
            it = tmp_it;
    }

    Iterator() = default;
    Iterator(const dynamic_pgm_type *p, const internal_iterator it) : super(p), it(it), initialized(), queue() {};

public:

    using difference_type = ssize_t;
    using value_type = const Item;
    using pointer = const Item *;
    using reference = const Item &;
    using iterator_category = std::forward_iterator_tag;

    Iterator &operator++() {
        lazy_initialize_queue();
        advance();
        return *this;
    }

    Iterator operator++(int) {
        Iterator i(it);
        ++i;
        return i;
    }

    reference operator*() const { return *it; }
    pointer operator->() const { return &*it; };
    bool operator==(const Iterator &rhs) const { return it == rhs.it; }
    bool operator!=(const Iterator &rhs) const { return it != rhs.it; }
};

#pragma pack(push, 1)

template<typename K, typename V, typename PGMType, uint8_t MinIndexedLevel>
class DynamicPGMIndex<K, V, PGMType, MinIndexedLevel>::ItemA {
    friend class DynamicPGMIndex;
    static V tombstone;

    ItemA() = default;
    explicit ItemA(const K &key) : first(key), second(tombstone) {}
    explicit ItemA(const K &key, const V &value) : first(key), second(value) {}

    bool deleted() const { return this->second == tombstone; }

public:
    K first;
    V second;

    operator K() const { return first; }
};

template<typename K, typename V, typename PGMType, uint8_t MinIndexedLevel>
V DynamicPGMIndex<K, V, PGMType, MinIndexedLevel>::ItemA::tombstone = new std::remove_pointer_t<V>();

template<typename K, typename V, typename PGMType, uint8_t MinIndexedLevel>
class DynamicPGMIndex<K, V, PGMType, MinIndexedLevel>::ItemB {
    friend class DynamicPGMIndex;
    bool flag;

    ItemB() = default;
    explicit ItemB(const K &key) : flag(true), first(key), second() {}
    explicit ItemB(const K &key, const V &value) : flag(false), first(key), second(value) {}

    bool deleted() const { return flag; }

public:
    K first;
    V second;

    operator K() const { return first; }
};

#pragma pack(pop)

}