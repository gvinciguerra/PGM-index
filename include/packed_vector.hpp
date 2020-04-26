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

#include <vector>
#include <stdexcept>
#include <immintrin.h>

/**
 *  A container for integers whose bit-length is fixed. The size is not necessarily the one of standard types
 *  (e.g. 8, 16, 32, 64 bits).
 */
class PackedVector {
    size_t n;
    uint8_t bpi;
    std::vector<uint64_t> data;
    static constexpr uint8_t w = 64;

    inline void write(size_t i, uint64_t x) {
        const uint64_t u = 1;
        x &= (u << bpi) - 1;
        auto lo = i * bpi;
        auto hi = (i + 1) * bpi - 1;
        if (lo / w == hi / w) {
            data[hi / w] &= ~(((u << (hi - lo + 1)) - 1) << (lo % w));
            data[hi / w] |= x << (lo % w);
        } else {
            data[lo / w] = (data[lo / w] & ((u << (lo % w)) - 1)) | (x << (lo % w));
            data[hi / w] = (data[hi / w] & ~((u << ((hi + 1) % w)) - 1)) | (x >> (w - (lo % w)));
        }
    }

public:
    
    /// Returns the length in bits of an integer x.
    static inline uint8_t bit_length(uint64_t x) {
        return x == 0 ? 1 : 64 - __builtin_clzl(x);
    }

    /**
     * Constructs an empty container.
     */
    PackedVector() : data(), n(0), bpi(0) {}

    /**
     * Constructs a new container with the specified bit-length from the specified vector.
     * @tparam T the type of the elements in the vector
     * @param in the vector to be used to initialize the elements of the container with
     * @param bpi the bit-length to use for the container
     */
    template<typename T>
    PackedVector(std::vector<T> &in, uint8_t bpi) : data(bpi * in.size() / w + 2), n(in.size()), bpi(bpi) {
        if (bpi < 1 || bpi > 63)
            throw std::invalid_argument("!(0 < bpi < 64)");
        for (size_t i = 0; i < n; ++i)
            write(i, in[i]);
    }

    /**
     * Constructs a new container with the contents of the range [@p first, @p last).
     * @tparam Iter the type of the iterator
     * @param first the iterator to the first element
     * @param last the iterator to the last element
     * @param bpi the bit-length to use for the container
     */
    template<typename Iter>
    PackedVector(Iter first, Iter last, uint8_t bpi) : data(bpi * (last - first) / w + 2), n(last - first), bpi(bpi) {
        if (bpi < 1 || bpi > 63)
            throw std::invalid_argument("!(0 < bpi < 64)");
        for (size_t i = 0; i < n; ++i, ++first)
            write(i, *first);
    }

    /**
     * Constructs a new container from the specified vector, fixing the bit-length to the one of the largest value.
     * @tparam T the type of the elements in the vector
     * @param in the vector to be used to initialize the elements of the container with
     */
    template<typename T>
    explicit PackedVector(std::vector<T> &in)
        : PackedVector(in, PackedVector::bit_length(*std::max_element(in.begin(), in.end()))) {}

    /**
     * Returns the specified element. No bounds checking is performed.
     * @param i position of the element to return
     * @return the requested element
     */
    inline uint64_t operator[](size_t i) const {
#ifdef __BMI2__
        auto lo = i * bpi;
        auto hi = (i + 1) * bpi - 1;
        auto word_i = lo / w;
        auto carryover_length = lo % w + bpi > w ? (hi + 1) % w : 0;
        auto bits = _bextr_u64(data[word_i], lo % w, bpi);
        auto carryover_bits = _bzhi_u64(data[word_i + 1], carryover_length);
        auto res = (carryover_bits << (w - (lo % w))) | bits;
        return res;
#else
        auto lo = i * bpi;
        auto hi = (i + 1) * bpi - 1;
        const uint64_t u = 1;
        if (lo / w == hi / w)
            return (data[hi / w] >> (lo % w)) & ((u << (hi - lo + 1)) - 1);
        else
            return (data[lo / w] >> (lo % w)) | (data[hi / w] & ((u << ((hi + 1) % w)) - 1)) << (w - (lo % w));
#endif
    }

    PackedVector &operator=(const PackedVector &v) = default;

    /**
     * Returns the number of elements in the container.
     * @return the number of elements in the container
     */
    size_t size() const {
        return n;
    }

    /**
     * Returns the size in bytes of the container.
     * @return the size in bytes of the container
     */
    size_t size_in_bytes() const {
        return bpi * size() / 8;
    }

};