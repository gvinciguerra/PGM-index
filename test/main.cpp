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

#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "pgm_index.hpp"
#include "pgm_index_compressed.hpp"
#include <type_traits>

TEMPLATE_TEST_CASE("Segmentation algorithm", "", float, uint32_t) {
    const auto error = 64;
    std::vector<TestType> data(1000000);
    std::random_device rnd_device;
    std::mt19937 engine{rnd_device()};

    if constexpr (std::is_floating_point<TestType>()) {
        std::lognormal_distribution<TestType> dist(0, 0.5);
        std::generate(data.begin(), data.end(), [&dist, &engine]() { return dist(engine); });
    }
    else
        std::generate(data.begin(), data.end(), std::rand);

    std::sort(data.begin(), data.end());
    auto segments = Segmentation<TestType, 0>::build_segments(data, error);

    for (auto i = 0; i < data.size(); ++i) {
        if (i == 0 || data[i] != data[i - 1]) {
            auto key = data[i];
            auto pos = std::lower_bound(segments.cbegin(), segments.cend(), key);
            auto segment = pos == segments.cbegin() || key == pos->key ? pos : pos - 1;
            double e = std::fabs(double(i) - std::min(data.size() - 1, segment->operator()(key)));
            REQUIRE(e <= error);
        }
    }
}

using S = Segmentation<uint32_t, 32>;
TEMPLATE_TEST_CASE("PGM-index", "", (BinarySearchStrategy<S>), (TreeStrategy<S, 32>), (RecursiveStrategy<S, 32>)) {
    std::vector<uint32_t> data(1000000);
    std::generate(data.begin(), data.end(), [] { return std::rand() % 10000; });
    std::sort(data.begin(), data.end());
    PGMIndex<uint32_t, 32, TestType> pgm_index(data);

    for (auto i = 1; i <= 1000; ++i) {
        auto q = data[std::rand() % data.size()];
        auto approx_range = pgm_index.find_approximate_position(q);
        auto lo = data.cbegin() + approx_range.lo;
        auto hi = data.cbegin() + approx_range.hi;
        auto k = std::lower_bound(lo, hi, q);
        REQUIRE(*k == q);
    }
}

TEST_CASE("Compressed PGM-index") {
    std::vector<uint32_t> data(1000000);
    std::generate(data.begin(), data.end(), [] { return std::rand() % 10000; });
    std::sort(data.begin(), data.end());
    CompressedPGMIndex<uint32_t, 32, 32> compressed_pgm_index(data);

    for (auto i = 1; i <= 1000; ++i) {
        auto q = data[std::rand() % data.size()];
        auto approx_range = compressed_pgm_index.find_approximate_position(q);
        auto lo = data.cbegin() + approx_range.lo;
        auto hi = data.cbegin() + approx_range.hi;
        auto k = std::lower_bound(lo, hi, q);
        REQUIRE(*k == q);
    }
}