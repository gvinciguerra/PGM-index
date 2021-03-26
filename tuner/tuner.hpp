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

#include "benchmark.hpp"
#include "pgm/pgm_index.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

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
        this->first_key = data[0];
        this->build(data.begin(), data.end(), epsilon, PGM_EPSILON_RECURSIVE, this->segments, this->levels_offsets);
    }

    pgm::ApproxPos search(const K &key) const {
        auto k = std::max(this->first_key, key);
        auto it = this->segment_for_key(k);
        auto pos = std::min<size_t>((*it)(k), std::next(it)->intercept);
        auto lo = PGM_SUB_EPS(pos, epsilon);
        auto hi = PGM_ADD_EPS(pos, epsilon, this->n);
        return {pos, lo, hi};
    }
};

/*------- INDEX STATS -------*/

struct IndexStats {
    size_t epsilon;
    size_t segments_count;
    size_t bytes;
    size_t lookup_ns;
    size_t lookup_ns_std;
    size_t construction_ns;

    IndexStats() = default;

    template<typename K>
    IndexStats(const std::vector<K> &data, const std::vector<K> &queries, size_t epsilon) : epsilon(epsilon) {
        auto start = timer::now();
        MockPGMIndex<K> pgm(data, epsilon);
        auto end = timer::now();
        construction_ns = size_t(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        double avg_time = 0;
        double var_time = 0;
        const int repetitions = 5;

        for (int repetition = 1; repetition < repetitions; ++repetition) {
            auto t0 = timer::now();

            uint64_t cnt = 0;
            for (auto &q : queries) {
                auto range = pgm.search(q);
                auto lo = data.begin() + range.lo;
                auto hi = data.begin() + range.hi;
                cnt += std::distance(data.begin(), std::lower_bound(lo, hi, q));
            }
            [[maybe_unused]] volatile auto tmp = cnt;

            auto t1 = timer::now();
            auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() / queries.size();
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
    std::printf("%-19zu %-19.2f %-19.2f %-19s", stats.epsilon, stats.construction_ns * 1.e-9, kib, query_time.c_str());
    if (verbose) {
        auto bounds = "(" + std::to_string(lo_eps) + ", " + std::to_string(hi_eps) + ")";
        std::printf("\t↝ search space=%-15s", bounds.c_str());
    }
    std::printf("\n");
}

template<typename K>
void minimize_space_given_time(size_t max_time, double tolerance, const std::vector<K> &data,
                               const std::vector<K> &queries, size_t lo_eps, size_t hi_eps, bool verbose) {
    auto latency = 82.1;
    auto cache_line = cache_line_size();
    auto block_size = cache_line / sizeof(K);
    auto eps_start = std::clamp(size_t(block_size * std::pow(2., max_time / latency - 1.)), lo_eps, hi_eps);

    const size_t starting_i = 2048;
    auto i = starting_i;
    auto lo = lo_eps;
    auto hi = hi_eps;

    std::vector<IndexStats> all_stats;
    all_stats.emplace_back(data, queries, eps_start);
    minimize_time_logging(all_stats.back(), verbose, eps_start, eps_start);

    if (all_stats.back().lookup_ns < max_time) {
        while (eps_start + (i << 1) < hi_eps && all_stats.back().lookup_ns < max_time * (1 + tolerance)) {
            i <<= 1;
            all_stats.emplace_back(data, queries, eps_start + i);
            lo = eps_start + i / 2;
            hi = eps_start + i;
            minimize_time_logging(all_stats.back(), verbose, lo, hi);
        }
        lo = i <= (starting_i << 1) ? eps_start : eps_start + i / 2;
    } else {
        while (eps_start > (i << 1) + lo_eps && all_stats.back().lookup_ns > max_time * (1 - tolerance)) {
            i <<= 1;
            all_stats.emplace_back(data, queries, eps_start - i);
            lo = eps_start - i;
            hi = eps_start - i / 2;
            minimize_time_logging(all_stats.back(), verbose, lo, hi);
        }
        hi = i <= (starting_i << 1) ? eps_start : eps_start - i / 2;
    }

    while (hi - lo > cache_line / 2) {
        i = (hi + lo) / 2;
        all_stats.emplace_back(data, queries, i);
        if (all_stats.back().lookup_ns > max_time)
            hi = i;
        else
            lo = i + 1;
        minimize_time_logging(all_stats.back(), verbose, lo, hi);
    }

    auto pred = [&](const IndexStats &a) { return a.lookup_ns <= max_time * (1 + tolerance); };
    if (!std::any_of(all_stats.cbegin(), all_stats.cend(), pred)) {
        std::printf("It is not possible to satisfy the given constraint. Increase the maximum time.");
        std::exit(1);
    }

    auto cmp = [&](const IndexStats &a, const IndexStats &b) { return !pred(b) || (pred(a) && a.bytes < b.bytes); };
    auto best = std::min_element(all_stats.cbegin(), all_stats.cend(), cmp);
    std::printf("%s\n", std::string(80, '-').c_str());
    std::printf("Set epsilon to %zu for an index of %zu bytes\n", best->epsilon, best->bytes);
}

template<typename K>
void minimize_time_given_space(size_t max_space, double tolerance, const std::vector<K> &data,
                               const std::vector<K> &queries, size_t lo_eps, size_t hi_eps, bool verbose) {
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

        all_stats.emplace_back(data, queries, mid);
        auto &stats = all_stats.back();
        auto kib = stats.bytes / double(1u << 10u);
        auto query_time = std::to_string(stats.lookup_ns) + "±" + std::to_string(stats.lookup_ns_std);
        std::printf("%-19zu %-19.2f %-19.2f %-19s", mid, stats.construction_ns * 1.e-9, kib, query_time.c_str());
        if (verbose) {
            auto s = "(" + std::to_string(lo) + ", " + std::to_string(hi) + ")";
            std::printf("\t↝ search space=%-15s \ts(ε)=%.0fε^%.2f \tε guess=%zu", s.c_str(), a, b, guess);
        }
        std::printf("\n");
        std::fflush(stdout);

        if (stats.bytes <= max_space)
            hi = mid;
        else
            lo = mid + 1;
    } while (lo < hi && std::abs(all_stats.back().bytes - (double) max_space) > max_space * tolerance);

    std::printf("%s\n", std::string(80, '-').c_str());
    std::printf("Set epsilon to %zu\n", lo);
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

size_t cache_line_size() {
    FILE *p = std::fopen("/sys/devices/system/cpu/cpu0/cache/index0/coherency_line_size", "r");
    unsigned int i = 0;
    if (p) {
        std::fscanf(p, "%u", &i);
        std::fclose(p);
    }
    return i;
}

#else
#error Unrecognized platform
#endif