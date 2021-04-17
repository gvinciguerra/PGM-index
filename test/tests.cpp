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

#include "catch.hpp"
#include "pgm/morton_nd.hpp"
#include "pgm/pgm_index.hpp"
#include "pgm/pgm_index_dynamic.hpp"
#include "pgm/pgm_index_variants.hpp"
#include "pgm/piecewise_linear_model.hpp"
#include "utils.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

template <typename Index, typename Data>
void test_index(const Index &index, const Data &data) {
    auto rand = std::bind(std::uniform_int_distribution<size_t>(0, data.size() - 1), std::mt19937{42});

    for (auto i = 1; i <= 10000; ++i) {
        auto q = data[rand()];
        auto range = index.search(q);
        auto lo = data.begin() + range.lo;
        auto hi = data.begin() + range.hi;
        auto k = std::lower_bound(lo, hi, q);
        REQUIRE(*k == q);
    }

    // Test elements outside range
    auto q = data.back() + 42;
    auto range = index.search(q);
    auto lo = data.begin() + range.lo;
    auto hi = data.begin() + range.hi;
    REQUIRE(std::lower_bound(lo, hi, q) == data.end());

    q = 0;
    range = index.search(q);
    lo = data.begin() + range.lo;
    hi = data.begin() + range.hi;
    REQUIRE(std::lower_bound(lo, hi, q) == data.begin());
}

TEMPLATE_TEST_CASE("Segmentation algorithm", "", float, double, uint32_t, uint64_t) {
    auto epsilon = GENERATE(32, 64, 128);
    auto data = generate_data<TestType>(1000000);
    auto segments = pgm::internal::make_segmentation(data.begin(), data.end(), epsilon);
    auto it = segments.begin();
    auto [slope, intercept] = it->get_floating_point_segment(it->get_first_x());

    for (auto i = 0u; i < data.size(); ++i) {
        if (i != 0 && data[i] == data[i - 1])
            continue;
        if (std::next(it) != segments.end() && std::next(it)->get_first_x() <= data[i]) {
            ++it;
            std::tie(slope, intercept) = it->get_floating_point_segment(it->get_first_x());
        }

        auto pos = (data[i] - it->get_first_x()) * slope + intercept;
        auto e = std::fabs(i - pos);
        REQUIRE(e <= epsilon + 1);
    }
}

TEMPLATE_TEST_CASE_SIG("PGM-index", "",
                       ((typename T, size_t E1, size_t E2), T, E1, E2),
                       (uint32_t, 16, 0), (uint32_t, 32, 0), (uint32_t, 64, 0),
                       (uint64_t, 16, 4), (uint64_t, 32, 4), (uint64_t, 64, 4),
                       (uint64_t, 4, 16), (uint64_t, 4, 32), (uint64_t, 4, 64)) {
    auto data = generate_data<T>(3000000);
    pgm::PGMIndex<T, E1, E2> index(data.begin(), data.end());
    test_index(index, data);
}

TEMPLATE_TEST_CASE_SIG("Compressed PGM-index", "", ((size_t E), E), 8, 32, 128) {
    auto data = generate_data<uint32_t>(3000000);
    pgm::CompressedPGMIndex<uint32_t, E> index(data);
    test_index(index, data);
}

TEMPLATE_TEST_CASE_SIG("Bucketing PGM-index", "", ((size_t E), E), 8, 32, 128) {
    auto data = generate_data<uint32_t>(3000000);
    auto top_level_size = GENERATE(256, 1024, 4096);
    pgm::BucketingPGMIndex<uint32_t, E> index(data.begin(), data.end(), top_level_size);
    test_index(index, data);
}


TEST_CASE("Bucketing PGM-index edge case", "") {
    std::vector<uint32_t> data(3000000);
    auto top_level_size = GENERATE(256, 1024, 4096);
    pgm::BucketingPGMIndex<uint32_t, 16> index(data.begin(), data.end(), top_level_size);
    test_index(index, data);
}

TEMPLATE_TEST_CASE_SIG("Elias-Fano PGM-index", "", ((size_t E), E), 8, 32, 128) {
    auto data = generate_data<uint32_t>(3000000);
    pgm::EliasFanoPGMIndex<uint32_t, E> index(data.begin(), data.end());
    test_index(index, data);
}

TEMPLATE_TEST_CASE_SIG("Mapped PGM-index", "", ((size_t E), E), 8, 32, 128) {
    std::string tmp_filename = "tmp.mapped.pgm";
    auto data = generate_data<uint32_t>(500000);
    auto random_query = std::bind(std::uniform_int_distribution<uint32_t>(data.front(), data.back()), std::mt19937{42});

    {
        pgm::MappedPGMIndex<uint32_t, E> index(data.begin(), data.end(), tmp_filename);
        for (auto i = 1; i <= 5000; ++i) {
            auto q = random_query();
            auto lb = std::lower_bound(data.begin(), data.end(), q);
            auto ub = std::upper_bound(data.begin(), data.end(), q);
            REQUIRE(std::distance(index.begin(), index.lower_bound(q)) == std::distance(data.begin(), lb));
            REQUIRE(std::distance(index.begin(), index.upper_bound(q)) == std::distance(data.begin(), ub));
        }
    }

    {
        pgm::MappedPGMIndex<uint32_t, E> index(tmp_filename);
        for (auto i = 1; i <= 5000; ++i) {
            auto q = random_query();
            REQUIRE(index.count(q) == (size_t) std::count(data.begin(), data.end(), q));
        }
    }

    std::remove(tmp_filename.c_str());
}

TEMPLATE_TEST_CASE_SIG("Dynamic PGM-index", "",
                       ((typename V, uint8_t MinIndexedLevel), V, MinIndexedLevel),
                       (uint32_t*, 8), (uint32_t, 10), (uint32_t*, 16), (uint32_t, 20)) {
    V time = 0;
    auto rand = std::bind(std::uniform_int_distribution<uint32_t>(0, 1000000000), std::mt19937{42});
    auto gen = [&] { return std::pair<uint32_t, V>{rand(), ++time}; };

    std::vector<std::pair<uint32_t, V>> bulk(GENERATE(0, 10, 1000, 1000000));
    std::generate(bulk.begin(), bulk.end(), gen);
    std::sort(bulk.begin(), bulk.end());
    bulk.erase(std::unique(bulk.begin(), bulk.end(), [](auto &a, auto &b) { return a.first == b.first; }), bulk.end());

    using PGMType = pgm::PGMIndex<uint32_t>;
    pgm::DynamicPGMIndex<uint32_t, V, PGMType, MinIndexedLevel> pgm(bulk.begin(), bulk.end());
    std::map<uint32_t, V> map(bulk.begin(), bulk.end());

    // Test initial state
    auto it1 = pgm.begin();
    for (auto[k, v] : map) {
        REQUIRE(it1->first == k);
        REQUIRE(it1->second == v);
        ++it1;
    }

    // Test lower bound
    for (size_t i = 0; i < std::min<size_t>(1000, bulk.size()); ++i) {
        auto q = bulk[rand() % bulk.size()];
        auto c = pgm.count(q.first);
        auto it = pgm.lower_bound(q.first);
        REQUIRE(c == 1);
        REQUIRE(it->first == q.first);
    }

    // Overwrite some elements
    ++time;
    for (size_t i = 0; i < std::min<size_t>(10000, bulk.size()); ++i, ++time) {
        pgm.insert_or_assign(bulk[i].first, time);
        map.insert_or_assign(bulk[i].first, time);
    }

    // Insert new elements
    for (size_t i = 0; i < 10000; ++i) {
        auto[k, v] = gen();
        pgm.insert_or_assign(k, v);
        map.insert_or_assign(k, v);
    }
    REQUIRE(pgm.size() == map.size());

    // Test for most recent values
    for (size_t i = 0; i < std::min<size_t>(10000, bulk.size()); ++i) {
        auto q = bulk[i];
        auto it = pgm.lower_bound(q.first);
        REQUIRE(it->first == q.first);
        REQUIRE(it->second > q.second);
        REQUIRE(it->second == map.lower_bound(q.first)->second);
    }

    // Delete some elements
    for (size_t i = 10; i < std::min<size_t>(500, bulk.size()); ++i) {
        pgm.erase(bulk[i].first);
        map.erase(bulk[i].first);
    }

    // Check if elements are deleted
    auto end = pgm.end();
    for (size_t i = 10; i < std::min<size_t>(500, bulk.size()); ++i) {
        auto it = pgm.find(bulk[i].first);
        REQUIRE(it == end);
    }
    REQUIRE(pgm.size() == map.size());

    // Test iterator
    auto it = pgm.begin();
    for (auto[k, v] : map) {
        REQUIRE(it->first == k);
        REQUIRE(it->second == v);
        ++it;
    }
}

#ifdef MORTON_ND_BMI2_ENABLED

TEMPLATE_TEST_CASE_SIG("Multidimensional PGM-index", "",
                       ((typename T, uint8_t D), T, D),
                       (uint32_t, 2), (uint32_t, 3), (uint64_t, 2), (uint64_t, 3), (uint64_t, 4)) {
    auto u = 1ull << (std::numeric_limits<T>::digits / D - 2);
    auto rand = std::bind(std::uniform_int_distribution<T>(0, u), std::mt19937{42});
    auto rand_tuple = [&] { return make_rand_tuple(rand, std::make_index_sequence<D>()); };

    std::vector<decltype(rand_tuple())> data(1000000);
    std::generate(data.begin(), data.end(), rand_tuple);
    pgm::MultidimensionalPGMIndex<D, T, 16> pgm(data.begin(), data.end());

    for (auto p: data)
        REQUIRE(pgm.contains(p));

    for (int i = 0; i < 500; ++i) {
        auto min = rand_tuple();
        auto max = min + rand_tuple();
        auto count = std::distance(pgm.range(min, max), pgm.end());
        auto expected_count = 0;
        for (auto &x : data)
            expected_count += box_contains(min, max, x);
        REQUIRE(count == expected_count);
    }
}

#endif