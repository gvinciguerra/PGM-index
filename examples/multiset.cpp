/*
 * This example shows how to use the PGM-index to implement classical query operations on a sorted multiset.
 * Compile with:
 *   g++ multiset.cpp -std=c++17 -I../include -o multiset
 * Run with:
 *   ./multiset
 */

#include <vector>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <algorithm>
#include "pgm/pgm_index.hpp"

template<typename K>
class PGMMultiset {
    std::vector<K> data;
    pgm::PGMIndex<K, 32, 4, float> pgm;

public:

    explicit PGMMultiset(const std::vector<K> &data) : data(data), pgm(data.begin(), data.end()) {}

    bool contains(const K x) const {
        auto range = pgm.search(x);
        return std::binary_search(data.begin() + range.lo, data.begin() + range.hi, x);
    }

    auto lower_bound(const K x) const {
        auto range = pgm.search(x);
        return std::lower_bound(data.begin() + range.lo, data.begin() + range.hi, x);
    }

    auto upper_bound(const K x) const {
        auto range = pgm.search(x);
        auto it = std::upper_bound(data.begin() + range.lo, data.begin() + range.hi, x);
        auto step = 1ull;
        while (it + step < end() && *(it + step) == x)  // exponential search to skip duplicates
            step *= 2;
        return std::upper_bound(it + (step / 2), std::min(it + step, end()), x);
    }

    size_t count(const K x) const {
        auto lb = lower_bound(x);
        if (lb == end() || *lb != x)
            return 0;
        return std::distance(lb, upper_bound(x));
    }

    auto begin() const { return data.cbegin(); }
    auto end() const { return data.cend(); }
};

template<typename T, typename Iter>
void print_iterator(Iter first, Iter last) {
    std::cout << "[";
    if (first != last) {
        std::copy(first, std::prev(last), std::ostream_iterator<T>(std::cout, " "));
        std::cout << *std::prev(last);
    }
    std::cout << "]" << std::endl;
}

int main() {
    std::vector<int> data(100);
    std::generate(data.begin(), data.end(), []() { return std::rand() % 100; });
    std::sort(data.begin(), data.end());

    std::cout << "data = ";
    print_iterator<int>(data.begin(), data.end());

    PGMMultiset multiset(data);
    std::cout << "count(15) = " << multiset.count(15) << std::endl;
    std::cout << "contains(1) = " << (multiset.contains(1) ? "true" : "false") << std::endl;
    std::cout << "lower_bound(42) = " << *multiset.lower_bound(42) << std::endl;
    std::cout << "upper_bound(70) = " << *multiset.upper_bound(70) << std::endl;

    std::cout << "Range search [50, 60] = ";
    print_iterator<int>(multiset.lower_bound(50), multiset.upper_bound(60));

    return 0;
}