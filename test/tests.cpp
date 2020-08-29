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

#include <map>
#include <random>
#include <functional>
#include <type_traits>
#include "catch.hpp"
#include "pgm/pgm_index.hpp"
#include "pgm/pgm_index_dynamic.hpp"
#include "pgm/pgm_index_compressed.hpp"

TEMPLATE_TEST_CASE("Segmentation algorithm", "", float, double, uint32_t, uint64_t) {
    const auto epsilon = GENERATE(32, 64, 128);
    std::vector<TestType> data(1000000);
    std::mt19937 engine(42);
    using RandomFunction = std::function<TestType()>;

    if constexpr (std::is_floating_point<TestType>()) {
        RandomFunction lognormal = std::bind(std::lognormal_distribution<TestType>(0, 0.5), engine);
        RandomFunction exponential = std::bind(std::exponential_distribution<TestType>(1.2), engine);
        auto rand = GENERATE_COPY(as<RandomFunction>{}, lognormal, exponential);
        std::generate(data.begin(), data.end(), rand);
    } else {
        RandomFunction uniform_dense = std::bind(std::uniform_int_distribution<TestType>(0, 10000), engine);
        RandomFunction uniform_sparse = std::bind(std::uniform_int_distribution<TestType>(0, 10000000), engine);
        RandomFunction binomial = std::bind(std::binomial_distribution<TestType>(50000), engine);
        RandomFunction geometric = std::bind(std::geometric_distribution<TestType>(0.8), engine);
        auto rand = GENERATE_COPY(as<RandomFunction>{}, uniform_dense, uniform_sparse, binomial, geometric);
        std::generate(data.begin(), data.end(), rand);
    }

    std::sort(data.begin(), data.end());
    auto segments = pgm::internal::make_segmentation(data.begin(), data.end(), epsilon);
    auto it = segments.begin();
    auto [slope, intercept] = it->get_floating_point_segment(it->get_first_x());

    for (auto i = 0; i < data.size(); ++i) {
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
    std::vector<T> data(2000000);
    std::mt19937 engine(42);

    using RandomFunction = std::function<T()>;
    RandomFunction uniform_dense = std::bind(std::uniform_int_distribution<T>(0, 10000), engine);
    RandomFunction uniform_sparse = std::bind(std::uniform_int_distribution<T>(0, 10000000), engine);
    RandomFunction binomial = std::bind(std::binomial_distribution<T>(50000), engine);
    RandomFunction geometric = std::bind(std::geometric_distribution<T>(0.8), engine);
    auto rand = GENERATE_COPY(as<RandomFunction>{}, uniform_dense, uniform_sparse, binomial, geometric);

    std::generate(data.begin(), data.end(), rand);
    std::sort(data.begin(), data.end());
    pgm::PGMIndex<T, E1, E2> pgm_index(data);

    for (auto i = 1; i <= 10000; ++i) {
        auto q = data[std::rand() % data.size()];
        auto range = pgm_index.search(q);
        auto lo = data.begin() + range.lo;
        auto hi = data.begin() + range.hi;
        auto k = std::lower_bound(lo, hi, q);
        REQUIRE(*k == q);
    }

    // Test elements outside range
    auto q = data.back() + 42;
    auto range = pgm_index.search(q);
    auto lo = data.begin() + range.lo;
    auto hi = data.begin() + range.hi;
    REQUIRE(std::lower_bound(lo, hi, q) == data.end());

    q = 0;
    range = pgm_index.search(q);
    lo = data.begin() + range.lo;
    hi = data.begin() + range.hi;
    REQUIRE(std::lower_bound(lo, hi, q) == data.begin());
}

TEST_CASE("Compressed PGM-index") {
    std::srand(42);
    std::vector<uint32_t> data(1000000);
    std::generate(data.begin(), data.end(), [] { return std::rand() % 10000; });
    std::sort(data.begin(), data.end());
    pgm::CompressedPGMIndex<uint32_t, 32, 32> compressed_pgm_index(data);

    for (auto i = 1; i <= 1000; ++i) {
        auto q = data[std::rand() % data.size()];
        auto range = compressed_pgm_index.search(q);
        auto lo = data.begin() + range.lo;
        auto hi = data.begin() + range.hi;
        auto k = std::lower_bound(lo, hi, q);
        REQUIRE(*k == q);
    }
}

TEMPLATE_TEST_CASE_SIG("Dynamic PGM-index", "",
                       ((typename V, uint8_t MinIndexedLevel), V, MinIndexedLevel),
                       (uint32_t*, 8), (uint32_t, 10), (uint32_t*, 16), (uint32_t, 20)) {
    V time = 0;
    std::srand(42);
    auto gen = [&time] { return std::pair<uint32_t, V>{std::rand() % 1000000000, ++time}; };

    std::vector<std::pair<uint32_t, V>> bulk(GENERATE(0, 10, 1000, 1000000));
    std::generate(bulk.begin(), bulk.end(), gen);
    std::sort(bulk.begin(), bulk.end());

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
        auto q = bulk[std::rand() % bulk.size()];
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