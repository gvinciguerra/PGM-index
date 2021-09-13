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

    q = std::numeric_limits<typename Data::value_type>::min();
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
                       (uint32_t, 8, 0), (uint32_t, 32, 0), (uint32_t, 128, 0),
                       (uint64_t, 8, 4), (uint64_t, 32, 4), (uint64_t, 128, 4),
                       (uint64_t, 256, 256), (uint64_t, 512, 512)) {
    auto data = generate_data<T>(2000000);
    pgm::PGMIndex<T, E1, E2> index(data.begin(), data.end());
    test_index(index, data);
}

TEMPLATE_TEST_CASE_SIG("Compressed PGM-index", "", ((size_t E), E), 8, 32, 128) {
    auto data = generate_data<uint32_t>(2000000);
    pgm::CompressedPGMIndex<uint32_t, E> index(data);
    test_index(index, data);
}

TEMPLATE_TEST_CASE_SIG("Bucketing PGM-index", "",
                       ((size_t E, size_t S), E, S), (4, 128), (8, 100), (4, 512), (8, 550)) {
    auto data = generate_data<uint32_t>(2000000);
    pgm::BucketingPGMIndex<uint32_t, E, S> index(data.begin(), data.end());
    test_index(index, data);
}

TEMPLATE_TEST_CASE_SIG("Bucketing PGM-index edge case", "",
                       ((size_t E, size_t S), E, S), (4, 128), (8, 100), (4, 512), (8, 550)) {
    std::vector<uint32_t> data(2000000);
    pgm::BucketingPGMIndex<uint32_t, E, S> index(data.begin(), data.end());
    test_index(index, data);
}

TEMPLATE_TEST_CASE_SIG("Elias-Fano PGM-index", "", ((size_t E), E), 8, 32, 128) {
    auto data = generate_data<uint32_t>(2000000);
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

TEMPLATE_TEST_CASE("Dynamic PGM-index", "", uint32_t*, uint32_t, std::string) {
    using time_type = uint32_t;
    auto make_key = std::bind(std::uniform_int_distribution<uint32_t>(0, 1000000000), std::mt19937{42});
    auto make_value = [&] {
        static time_type time = 0;
        if constexpr (std::is_same_v<TestType, std::string>) return std::to_string(++time);
        else return reinterpret_cast<TestType>(++time);
    };
    auto get_value = [](auto x) {
        if constexpr (std::is_same_v<TestType, std::string>) return (time_type) std::stoll(x);
        else return x;
    };
    auto gen = [&] { return std::pair<uint32_t, TestType>{make_key(), make_value()}; };

    std::vector<std::pair<uint32_t, TestType>> bulk(GENERATE(0, 10, 1000, 100000));
    std::generate(bulk.begin(), bulk.end(), gen);
    std::sort(bulk.begin(), bulk.end());
    bulk.erase(std::unique(bulk.begin(), bulk.end(), [](auto &a, auto &b) { return a.first == b.first; }), bulk.end());

    using PGMType = pgm::PGMIndex<uint32_t>;
    pgm::DynamicPGMIndex<uint32_t, TestType, PGMType> pgm(bulk.begin(), bulk.end(), GENERATE(2, 4, 8));
    std::map<uint32_t, TestType> map(bulk.begin(), bulk.end());

    // Test initial state
    auto it1 = pgm.begin();
    for (auto[k, v] : map) {
        REQUIRE(it1->first == k);
        REQUIRE(it1->second == v);
        ++it1;
    }

    // Test lower bound
    for (size_t i = 0; i < std::min<size_t>(1000, bulk.size()); ++i) {
        auto q = bulk[make_key() % bulk.size()];
        auto c = pgm.count(q.first);
        auto it = pgm.lower_bound(q.first);
        REQUIRE(c == 1);
        REQUIRE(it->first == q.first);
    }

    // Overwrite some elements
    for (size_t i = 0; i < std::min<size_t>(10000, bulk.size()); ++i) {
        auto v = make_value();
        pgm.insert_or_assign(bulk[i].first, v);
        map.insert_or_assign(bulk[i].first, v);
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
        REQUIRE(get_value(it->second) > get_value(q.second));
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

    // Test range
    for (int i = 0; i < 10; ++i) {
        auto lo = make_key();
        auto hi = lo + make_key() / 2;
        auto range_result = pgm.range(lo, hi);
        auto map_it = map.lower_bound(lo);
        for (auto[k, v] : range_result) {
            REQUIRE(k == map_it->first);
            REQUIRE(v == map_it->second);
            ++map_it;
        }
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