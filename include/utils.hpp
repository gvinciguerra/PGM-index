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

#include <cstddef>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <cstring>

std::vector<int64_t> read_data_csv(const std::string &filename,
                                   size_t max_lines = std::numeric_limits<size_t>::max()) {
    std::fstream in(filename);
    std::string line;
    std::vector<int64_t> data;
    data.reserve(max_lines == std::numeric_limits<size_t>::max() ? 1024 : max_lines);
    for (size_t i = 0; i < max_lines && std::getline(in, line); ++i) {
        int value;
        std::stringstream stringstream(line);
        stringstream >> value;
        data.push_back(value);
    }
    return data;
}

template<typename TypeIn, typename TypeOut>
std::vector<TypeOut> read_data_binary(const std::string &filename,
                                      size_t max_size = std::numeric_limits<size_t>::max()) {
    try {
        std::fstream in(filename, std::ios::in | std::ios::binary | std::ios::ate);
        in.exceptions(std::ios::failbit | std::ios::badbit);

        auto size = std::min(max_size, static_cast<size_t>(in.tellg() / sizeof(TypeIn)));
        std::vector<TypeIn> data(size);
        in.seekg(0);
        in.read((char *) data.data(), size * sizeof(TypeIn));
        if constexpr (std::is_same<TypeIn, TypeOut>::value)
            return data;

        return std::vector<TypeOut>(data.begin(), data.end());
    }
    catch (std::ios_base::failure &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << std::strerror(errno) << std::endl;
        exit(1);
    }
}

template<class RandomAccessIterator, class T>
bool interpolation_search(const RandomAccessIterator first, const RandomAccessIterator last, const T &value) {
    RandomAccessIterator lo = first;
    RandomAccessIterator hi = last - 1;

    if (value < *lo || value > *hi)
        return false;

    while (lo <= hi && value >= *lo && value <= *hi) {
        auto middle = lo + ((double) (hi - lo) / (*hi - *lo)) * (value - *lo);
        if (*middle < value)
            lo = middle + 1;
        else if (value < *middle)
            hi = middle - 1;
        else
            return true;
    }
    return false;
}

template<typename RandomAccessIterator, typename Compare, typename T>
RandomAccessIterator exponential_search(const RandomAccessIterator mid, const RandomAccessIterator first,
                                        const RandomAccessIterator last, const T &value, Compare comp) {
    auto lo = first;
    auto hi = last;
    auto i = 1;

    if (comp(*mid, value)) {
        while (mid + i < last && comp(*(mid + i), value))
            i <<= 1;
        return std::lower_bound(mid + (i / 2), std::min(mid + i, last), value, comp);
    } else {
        while (mid - i > first && !comp(*(mid - i), value))
            i <<= 1;
        return std::lower_bound(std::max(mid - i, first), mid - (i / 2), value, comp);
    }
}

/**
 * Returns the cache line size in bytes.
 * @return the cache line size in bytes
 */
size_t x86_cache_line() {
#define cpuid(func, ax, bx, cx, dx) __asm__ __volatile__ ("cpuid" : "=a"(ax), "=b"(bx), "=c"(cx), "=d"(dx) : "a"(func))
    size_t _, c;
    cpuid(0x80000006, _, _, c, _);
    return c & 0xFF;
}

template<class T>
void do_not_optimize(T const &value) {
    asm volatile("" : : "r,m"(value) : "memory");
}