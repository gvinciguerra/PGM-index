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

#include <cmath>
#include <cstdint>
#include <vector>
#include <chrono>
#include <numeric>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "pgm/pgm_index.hpp"

#define PGM_IGNORED_PARAMETER 1
#define PGM_EPSILON_RECURSIVE 4

template<typename K>
class MockPGMIndex : public pgm::PGMIndex<K, PGM_IGNORED_PARAMETER, PGM_EPSILON_RECURSIVE, double> {
    size_t epsilon;

public:

    using segment_type = typename pgm::PGMIndex<K, PGM_IGNORED_PARAMETER, PGM_EPSILON_RECURSIVE, double>::Segment;

    MockPGMIndex() = default;

    MockPGMIndex(const std::vector<K> &data, size_t epsilon) : epsilon(epsilon) {
        this->n = data.size();
        this->build(data.begin(), data.end(), epsilon, PGM_EPSILON_RECURSIVE,
                    this->segments, this->levels_sizes, this->levels_offsets);
    }

    pgm::ApproxPos search(const K &key) const {
        auto k = std::max(this->first_key, key);
        auto it = this->segment_for_key(k);
        auto pos = std::min<size_t>((*it)(k), std::next(it)->intercept);
        auto lo = PGM_SUB_EPS(pos, epsilon);
        auto hi = PGM_ADD_EPS(pos, epsilon, this->n);
        return {pos, lo, hi};
    }

    auto lower_bound(const std::vector<K> &data, K x) const {
        auto range = search(x);
        return std::lower_bound(data.begin() + range.lo, data.begin() + range.hi, x);
    }

};

/*------- INPUT FILE READING -------*/

std::vector<int64_t> read_data_csv(const std::string &file, size_t max_lines = std::numeric_limits<size_t>::max()) {
    std::fstream in(file);
    in.exceptions(std::ios::failbit | std::ios::badbit);
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
std::vector<TypeOut> read_data_binary(const std::string &file, size_t max_size = std::numeric_limits<size_t>::max()) {
    std::fstream in(file, std::ios::in | std::ios::binary | std::ios::ate);
    in.exceptions(std::ios::failbit | std::ios::badbit);

    auto size = std::min(max_size, static_cast<size_t>(in.tellg() / sizeof(TypeIn)));
    std::vector<TypeIn> data(size);
    in.seekg(0);
    in.read((char *) data.data(), size * sizeof(TypeIn));

    if constexpr (std::is_same<TypeIn, TypeOut>::value)
        return data;
    return std::vector<TypeOut>(data.begin(), data.end());
}

/*------- INDEX STATS -------*/

template<class T>
void do_not_optimize(T const &value) {
    asm volatile("" : : "r,m"(value) : "memory");
}

struct IndexStats {
    size_t epsilon;
    size_t segments_count;
    size_t bytes;
    size_t lookup_ns;
    size_t lookup_ns_std;
    size_t construction_ns;

    IndexStats() = default;

    template<typename K>
    IndexStats(const std::vector<K> &data, size_t epsilon) : epsilon(epsilon) {
        auto start = std::chrono::high_resolution_clock::now();
        MockPGMIndex<K> pgm(data, epsilon);
        auto end = std::chrono::high_resolution_clock::now();
        construction_ns = size_t(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        double avg_time = 0;
        double var_time = 0;
        const int repetitions = 5;
        const int queries = 100000;

        for (int repetition = 1; repetition < repetitions; ++repetition) {
            auto t0 = std::chrono::high_resolution_clock::now();

            for (int i = 1; i <= queries; ++i) {
                auto q = data[std::rand() % data.size()];
                do_not_optimize(pgm.lower_bound(data, q));
            }

            auto t1 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() / queries;
            auto next_avg_time = avg_time + (time - avg_time) / repetition;
            var_time += (time - avg_time) * (time - next_avg_time);
            avg_time = next_avg_time;
        }

        segments_count = pgm.segments_count();
        bytes = pgm.size_in_bytes();
        lookup_ns = size_t(avg_time);
        lookup_ns_std = size_t(std::sqrt(var_time / (repetitions - 1)));
    }
};

/*------- FUNCTION FITTING -------*/

/** Fits the coefficients (a,b) of a function f(ε)=aε^b. */
auto fit_segments_count_model(const std::vector<IndexStats> &all_index_stats) {
    auto count = 0;
    auto avg_y = 0.;
    auto avg_x = 0.;
    auto var_x = 0.;
    auto var_y = 0.;
    auto cov_xy = 0.;

    for (const auto &stats: all_index_stats) {
        count++;
        auto x = std::log(stats.epsilon);
        auto y = std::log(stats.segments_count);
        auto dx = x - avg_x;
        auto dy = y - avg_y;
        avg_x += dx / count;
        avg_y += dy / count;
        var_x += dx * (x - avg_x);
        var_y += dy * (y - avg_y);
        cov_xy += dx * (y - avg_y);
    }

    var_x /= (count - 1);
    cov_xy /= (count - 1);
    auto b = cov_xy / var_x;
    auto a = std::exp(avg_y - b * avg_x);
    return std::make_pair(a, b);
}

/*------- ROOT FINDING -------*/

double target_space_function(double epsilon, double a, double b, double max_space, double constants) {
    // s(ε) = c (m + (m-1) / (2ε-1), where m=aε^-b
    auto last_layer_segments = std::fmax(1., a * std::pow(epsilon, -b));
    auto total_segments = (2 * epsilon * last_layer_segments - 1) / (2 * epsilon - 1);
    return total_segments * constants - max_space;
}

double target_space_function_derivative(double epsilon, double a, double b, double, double constants) {
    // s'(ε) = (2c ε^(-b) (ab ε^b - a - 2abε)) / (2ε - 1)^2
    auto numerator = 2 * constants * std::pow(epsilon, -b) * (a * b + std::pow(epsilon, b) - a - 2 * a * b * epsilon);
    auto denominator = std::pow(2 * epsilon - 1, 2.);
    return numerator / denominator;
}

double guess_epsilon_space(double x, double a, double b, double max_space, double constants) {
    auto n = target_space_function(x, a, b, max_space, constants);
    auto d = target_space_function_derivative(x, a, b, max_space, constants);
    auto iterations = 0;

    while (std::abs(n / d) >= 0.0001 && iterations++ < 100) {
        n = target_space_function(x, a, b, max_space, constants);
        d = target_space_function_derivative(x, a, b, max_space, constants);
        x = x - n / d;
    }

    return std::fmax(0., x);
}

/*------- CORE FUNCTIONS -------*/

size_t cache_line_size();

void minimize_time_logging(const IndexStats &stats, bool verbose, size_t lo_eps, size_t hi_eps) {
    auto kib = stats.bytes / double(1u << 10u);
    auto query_time = std::to_string(stats.lookup_ns) + "±" + std::to_string(stats.lookup_ns_std);
    printf("%-19zu %-19.1f %-19.2f %-19s", stats.epsilon, stats.construction_ns * 1.e-9, kib, query_time.c_str());
    if (verbose) {
        auto bounds = "(" + std::to_string(lo_eps) + ", " + std::to_string(hi_eps) + ")";
        printf("\t↝ search space=%-15s", bounds.c_str());
    }
    printf("\n");
}

template<typename K>
void minimize_space_given_time(size_t max_time, double tolerance, const std::vector<K> &data,
                               size_t lo_eps, size_t hi_eps, bool verbose) {
    auto latency = 82.1;
    auto cache_line = cache_line_size();
    auto block_size = cache_line / sizeof(K);
    auto eps_start = std::clamp(size_t(block_size * std::pow(2., max_time / latency - 1.)), lo_eps, hi_eps);

    const size_t starting_i = 2048;
    auto i = starting_i;
    auto lo = lo_eps;
    auto hi = hi_eps;

    std::vector<IndexStats> all_stats;
    all_stats.emplace_back(data, eps_start);
    minimize_time_logging(all_stats.back(), verbose, eps_start, eps_start);

    if (all_stats.back().lookup_ns < max_time) {
        while (eps_start + (i << 1) < hi_eps && all_stats.back().lookup_ns < max_time * (1 + tolerance)) {
            i <<= 1;
            all_stats.emplace_back(data, eps_start + i);
            lo = eps_start + i / 2;
            hi = eps_start + i;
            minimize_time_logging(all_stats.back(), verbose, lo, hi);
        }
        lo = i <= (starting_i << 1) ? eps_start : eps_start + i / 2;
    } else {
        while (eps_start > (i << 1) + lo_eps && all_stats.back().lookup_ns > max_time * (1 - tolerance)) {
            i <<= 1;
            all_stats.emplace_back(data, eps_start - i);
            lo = eps_start - i;
            hi = eps_start - i / 2;
            minimize_time_logging(all_stats.back(), verbose, lo, hi);
        }
        hi = i <= (starting_i << 1) ? eps_start : eps_start - i / 2;
    }

    while (hi - lo > cache_line / 2) {
        i = (hi + lo) / 2;
        all_stats.emplace_back(data, i);
        if (all_stats.back().lookup_ns > max_time)
            hi = i;
        else
            lo = i + 1;
        minimize_time_logging(all_stats.back(), verbose, lo, hi);
    }

    auto pred = [&](const IndexStats &a) { return a.lookup_ns <= max_time * (1 + tolerance); };
    if (!std::any_of(all_stats.cbegin(), all_stats.cend(), pred)) {
        printf("It is not possible to satisfy the given constraint. Increase the maximum time.");
        exit(1);
    }

    auto cmp = [&](const IndexStats &a, const IndexStats &b) { return !pred(b) || (pred(a) && a.bytes < b.bytes); };
    auto best = std::min_element(all_stats.cbegin(), all_stats.cend(), cmp);
    auto time = std::accumulate(all_stats.cbegin(), all_stats.cend(), 0ul,
                                [](size_t result, const IndexStats &s) { return result + s.construction_ns; });

    printf("%s\n", std::string(80, '-').c_str());
    printf("%zu iterations for a total construction time of %.0f s\n", all_stats.size(), time * 1.e-9);
    printf("Set epsilon to %zu for an index of %zu bytes\n", best->epsilon, best->bytes);
}

template<typename K>
void minimize_time_given_space(size_t max_space, double tolerance, const std::vector<K> &data,
                               size_t lo_eps, size_t hi_eps, bool verbose) {
    const auto guess_steps_threshold = size_t(2 * std::log2(std::log2(hi_eps - lo_eps)));
    size_t guess_steps = 0;
    std::vector<IndexStats> all_stats;
    auto a = double(data.size() / 2);
    auto b = -1.;
    auto lo = lo_eps;
    auto hi = hi_eps;

    do {
        size_t guess = 0;
        size_t mid = (lo + hi) / 2;

        if (all_stats.size() >= 4 && guess_steps < guess_steps_threshold)
            std::tie(a,b) = fit_segments_count_model(all_stats);

        if (guess_steps < guess_steps_threshold) {
            auto constants = sizeof(typename MockPGMIndex<K>::segment_type);

            guess = size_t(guess_epsilon_space(100, a, -b, max_space, constants));
            guess = std::clamp(guess, lo + 1, hi - 1);

            auto bias_weight = guess_steps <= 1 ? 0 : double(guess_steps) / guess_steps_threshold;
            auto biased_guess = mid * bias_weight + guess * (1 - bias_weight);

            mid = size_t(biased_guess);
            guess_steps++;
        }

        all_stats.emplace_back(data, mid);
        auto &stats = all_stats.back();
        auto kib = stats.bytes / double(1u << 10u);
        auto query_time = std::to_string(stats.lookup_ns) + "±" + std::to_string(stats.lookup_ns_std);
        printf("%-19zu %-19.1f %-19.2f %-19s", mid, stats.construction_ns * 1.e-9, kib, query_time.c_str());
        if (verbose) {
            auto s = "(" + std::to_string(lo) + ", " + std::to_string(hi) + ")";
            printf("\t↝ search space=%-15s \ts(ε)=%.0fε^%.2f \tε guess=%zu", s.c_str(), a, b, guess);
        }
        printf("\n");
        fflush(stdout);

        if (stats.bytes <= max_space)
            hi = mid;
        else
            lo = mid + 1;
    } while (lo < hi && std::abs(all_stats.back().bytes - (double) max_space) > max_space * tolerance);

    auto time = std::accumulate(all_stats.cbegin(), all_stats.cend(), 0ull,
                                [](size_t r, const IndexStats &s) { return r + s.construction_ns; });
    printf("%s\n", std::string(80, '-').c_str());
    printf("Total construction time %.0f s\n", time * 1.e-9);
    printf("Set epsilon to %zu\n", lo);
}

/*------- cache_line_size() implementation (credits: https://stackoverflow.com/a/4049562) -------*/

#if defined(__APPLE__)

#include <sys/sysctl.h>
size_t cache_line_size() {
    size_t line_size = 0;
    size_t sizeof_line_size = sizeof(line_size);
    sysctlbyname("hw.cachelinesize", &line_size, &sizeof_line_size, 0, 0);
    return line_size;
}

#elif defined(_WIN32)

#include <stdlib.h>
#include <windows.h>
size_t cache_line_size() {
    size_t line_size = 0;
    DWORD buffer_size = 0;
    DWORD i = 0;
    SYSTEM_LOGICAL_PROCESSOR_INFORMATION * buffer = 0;

    GetLogicalProcessorInformation(0, &buffer_size);
    buffer = (SYSTEM_LOGICAL_PROCESSOR_INFORMATION *)malloc(buffer_size);
    GetLogicalProcessorInformation(&buffer[0], &buffer_size);

    for (i = 0; i != buffer_size / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION); ++i) {
        if (buffer[i].Relationship == RelationCache && buffer[i].Cache.Level == 1) {
            line_size = buffer[i].Cache.LineSize;
            break;
        }
    }

    free(buffer);
    return line_size;
}

#elif defined(linux)

#include <stdio.h>
size_t cache_line_size() {
    FILE *p = fopen("/sys/devices/system/cpu/cpu0/cache/index0/coherency_line_size", "r");
    unsigned int i = 0;
    if (p) {
        fscanf(p, "%u", &i);
        fclose(p);
    }
    return i;
}

#else
#error Unrecognized platform
#endif