// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2021 Giorgio Vinciguerra.
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

#include "catch.hpp"
#include <algorithm>
#include <cstddef>
#include <functional>
#include <iterator>
#include <random>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

template <typename T>
std::vector<T> generate_data(size_t n) {
    std::vector<T> data(n);
    std::mt19937 engine(42);

    using RandomFunction = std::function<T()>;
    if constexpr (std::is_floating_point<T>()) {
        RandomFunction lognormal = std::bind(std::lognormal_distribution<T>(0, 0.5), engine);
        RandomFunction exponential = std::bind(std::exponential_distribution<T>(1.2), engine);
        auto rand = GENERATE_COPY(as<RandomFunction>{}, lognormal, exponential);
        std::generate(data.begin(), data.end(), rand);
    } else {
        T min = 0;
        if constexpr (std::is_signed_v<T>)
            min = -10000;
        RandomFunction uniform_dense = std::bind(std::uniform_int_distribution<T>(min, 10000), engine);
        RandomFunction uniform_sparse = std::bind(std::uniform_int_distribution<T>(min, 10000000), engine);
        RandomFunction binomial = std::bind(std::binomial_distribution<T>(50000), engine);
        RandomFunction geometric = std::bind(std::geometric_distribution<T>(0.8), engine);
        auto rand = GENERATE_COPY(as<RandomFunction>{}, uniform_dense, uniform_sparse, binomial, geometric);
        std::generate(data.begin(), data.end(), rand);
    }

    std::sort(data.begin(), data.end());
    return data;
}

template<typename... T1, typename... T2, std::size_t... I>
constexpr auto add(const std::tuple<T1...> &t1, const std::tuple<T2...> &t2, std::index_sequence<I...>) {
    return std::tuple{std::get<I>(t1) + std::get<I>(t2)...}; // credits: https://stackoverflow.com/a/50815143
}

template<typename... T1, typename... T2>
constexpr auto operator+(const std::tuple<T1...>& t1, const std::tuple<T2...>& t2) {
    static_assert(sizeof...(T1) == sizeof...(T2));
    return add(t1, t2, std::make_index_sequence<sizeof...(T1)>{});
}

template<typename F, std::size_t... I>
auto make_rand_tuple(F f, std::index_sequence<I...>) {
    return std::make_tuple(f(I)...);
}

template<typename ... Ts, std::size_t ... Is>
bool box_contains_helper(const std::tuple<Ts...> &min,
                         const std::tuple<Ts...> &max,
                         const std::tuple<Ts...> &p,
                         const std::index_sequence<Is...> &) {
    return ((std::get<Is>(min) <= std::get<Is>(p) && std::get<Is>(p) <= std::get<Is>(max)) && ... );
}

template<typename ... Ts>
bool box_contains(const std::tuple<Ts...> &min, const std::tuple<Ts...> &max, const std::tuple<Ts...> &p) {
    return box_contains_helper(min, max, p, std::make_index_sequence<sizeof...(Ts)>{});
}

template<typename T1, typename T2>
constexpr bool box_contains(const std::pair<T1, T2> &min, const std::pair<T1, T2> &max, const std::pair<T1, T2> &p) {
    return box_contains(std::make_tuple(min.first, min.second),
                        std::make_tuple(max.first, max.second),
                        std::make_tuple(p.first, p.second));
}
