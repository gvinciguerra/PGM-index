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

#include "pgm_index.hpp"
#include "interpolation.h"

template<typename K>
class MockPGMIndex {
    size_t data_size;
    size_t error;

    struct SegmentData {
        double slope;
        double intercept;

        explicit SegmentData(const typename OptimalPiecewiseLinearModel<K, size_t>::CanonicalSegment &cs) {
            std::tie(slope, intercept) = cs.get_floating_point_segment(cs.get_first_x());
        }

        SegmentData(double slope, double intercept) : slope(slope), intercept(intercept) {};

        inline size_t operator()(const K &k) const {
            double pos = slope * k + intercept;
            return pos > 0. ? pos : 0ul;
        }
    };

    struct Layer {
        std::vector<K> segments_keys;
        std::vector<SegmentData> segments_data;

        inline size_t size() const {
            return segments_data.size();
        }

        template<typename S>
        explicit Layer(const S &segments) {
            segments_keys.reserve(segments.size());
            segments_data.reserve(segments.size());
            for (auto &cs : segments) {
                segments_keys.push_back(cs.get_first_x());
                segments_data.emplace_back(cs);
            }
        }
    };

    std::vector<Layer> layers;

public:

    MockPGMIndex(const std::vector<K> &data, size_t error) : data_size(data.size()), error(error) {
        std::list<Layer> tmp;
        tmp.emplace_front(make_segmentation(data.begin(), data.end(), error));

        while (tmp.front().size() > 1) {
            auto first = tmp.front().segments_keys.begin();
            auto last = tmp.front().segments_keys.end();
            tmp.emplace_front(make_segmentation(first, last, error));
        }

        layers = {std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end())};
    }

    inline K query(const std::vector<K> &data, K key) const {
        auto node_key = layers[0].segments_keys[0];
        size_t approx_pos = layers[0].segments_data[0](key - node_key);
        size_t pos = 0;

        for (auto it = layers.cbegin() + 1; it < layers.cend(); ++it) {
            auto layer_size = it->size();
            auto lo = SUB_ERR(approx_pos, error);
            auto hi = ADD_ERR(approx_pos, error + 1, layer_size);

            for (; lo <= hi && it->segments_keys[lo] <= key; ++lo);
            pos = lo - 1;

            node_key = it->segments_keys[pos];
            approx_pos = it->segments_data[pos](key - node_key);
            if (pos + 1 < it->size())
                approx_pos = std::min(approx_pos, (size_t) it->segments_data[pos + 1].intercept);
        }

        auto slope = layers.back().segments_data[pos].slope;
        auto intercept = layers.back().segments_data[pos].intercept;
        auto p = (size_t) std::fmax(0., slope * (key - node_key) + intercept);
        auto lo = SUB_ERR(p, error);
        auto hi = ADD_ERR(p, error, data_size);

        return *std::lower_bound(data.cbegin() + lo, data.cbegin() + hi + 1, key);
    }

    size_t size_in_bytes() const {
        auto total = 0;
        for (auto &l : layers)
            total += l.size();
        return total * (sizeof(SegmentData) + sizeof(K));
    }

    size_t segments_count() const {
        return layers.back().size();
    }

    size_t height() const {
        return layers.size();
    }

};

/*------- INDEX STATS -------*/

template<class T>
void do_not_optimize(T const &value) {
    asm volatile("" : : "r,m"(value) : "memory");
}

struct IndexStats {
    size_t error;
    size_t segments_count;
    size_t size_in_bytes;
    size_t lookup_time;
    size_t lookup_time_std;
    size_t construction_time;
    size_t height;

    IndexStats() {};

    template<typename K>
    IndexStats(std::vector<K> &data, size_t error) : error(error) {
        auto start = std::chrono::high_resolution_clock::now();
        MockPGMIndex<K> pgmIndex(data, error);
        auto end = std::chrono::high_resolution_clock::now();
        construction_time = size_t(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
        height = pgmIndex.height();

        double avg_time = 0;
        double var_time = 0;
        const int repetitions = 5;
        const int queries = 100000;

        for (int repetition = 1; repetition < repetitions; ++repetition) {
            auto t0 = std::chrono::high_resolution_clock::now();

            for (int i = 1; i <= queries; ++i) {
                auto q = data[std::rand() % data.size()];
                do_not_optimize(pgmIndex.query(data, q));
            }

            auto t1 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() / queries;

            double next_avg_time = avg_time + (time - avg_time) / repetition;
            var_time += (time - avg_time) * (time - next_avg_time);
            avg_time = next_avg_time;
        }

        segments_count = pgmIndex.segments_count();
        size_in_bytes = pgmIndex.size_in_bytes();
        lookup_time = size_t(avg_time);
        lookup_time_std = size_t(std::sqrt(var_time / (repetitions - 1)));
    }
};

/*------- FUNCTION FITTING -------*/

void segments_count_model(const alglib::real_1d_array &c, const alglib::real_1d_array &x, double &func, void *ptr) {
    func = c[0] * std::pow(x[0], -c[1]);
}

void segments_count_model_grad_c(const alglib::real_1d_array &c, const alglib::real_1d_array &x, double &func,
                                 alglib::real_1d_array &grad, void *ptr) {
    func = c[0] * std::pow(x[0], -c[1]);
    grad[0] = std::pow(x[0], -c[1]);
    grad[1] = -grad[0] * c[0] * std::log(x[0]);
}

auto fit_segments_count_model(const std::vector<IndexStats> &all_index_stats, alglib::real_1d_array c_space) {
    alglib::real_2d_array x;
    alglib::real_1d_array y_space;
    x.setlength(all_index_stats.size(), 1);
    y_space.setlength(all_index_stats.size());
    for (int r = 0; r < all_index_stats.size(); ++r) {
        x[r][0] = all_index_stats[r].error;
        y_space[r] = all_index_stats[r].segments_count;
    }

    alglib::real_1d_array s = "[1.0e+8, 1]";
    alglib::ae_int_t info;
    alglib::lsfitstate state;
    alglib::lsfitreport rep;
    alglib::lsfitcreatefg(x, y_space, c_space, true, state);
    alglib::lsfitsetcond(state, 0.00001, 500);
    alglib::lsfitsetscale(state, s);
    alglib::lsfitfit(state, segments_count_model, segments_count_model_grad_c);
    alglib::lsfitresults(state, info, c_space, rep);
    return std::make_pair(c_space, rep);
}

/*------- ROOT FINDING -------*/

double target_space_function(double error, double a, double b, double max_space, double constants) {
    // s(ε) = c (m + (m-1) / (2ε-1), where m=aε^-b
    double last_layer_segments = std::fmax(1., a * std::pow(error, -b));
    double total_segments = (2 * error * last_layer_segments - 1) / (2. * error - 1);
    return total_segments * constants - max_space;
}

double target_space_function_derivative(double error, double a, double b, double max_space, double constants) {
    // s'(ε) = (2c ε^(-b) (ab ε^b - a - 2abε)) / (2ε - 1)^2
    double numerator = 2 * constants * std::pow(error, -b) * (a * b + std::pow(error, b) - a - 2 * a * b * error);
    double denominator = std::pow(2 * error - 1, 2.);
    return numerator / denominator;
}

double guess_error_space(double x, double a, double b, double max_space, double constants) {
    double n = target_space_function(x, a, b, max_space, constants);
    double d = target_space_function_derivative(x, a, b, max_space, constants);

    while (std::abs(n / d) >= 0.0001) {
        n = target_space_function(x, a, b, max_space, constants);
        d = target_space_function_derivative(x, a, b, max_space, constants);
        x = x - n / d;
    }

    return std::fmax(0., x);
}

/*------- CORE FUNCTIONS -------*/

size_t x86_cache_line() {
    size_t _, c;
    __asm__ __volatile__("cpuid":"=a"(_), "=b"(_), "=c"(c), "=d"(_):"a"(0x80000006));
    return c & 0xFF;
}

void minimize_time_logging(IndexStats &stats, bool verbose, size_t lo_error, size_t hi_error) {
    auto kib = stats.size_in_bytes / double(1u << 10u);
    auto query_time = std::to_string(stats.lookup_time) + "±" + std::to_string(stats.lookup_time_std);
    printf("%-19zu %-19.1f %-19.2f %-19s", stats.error, stats.construction_time * 1.e-9, kib, query_time.c_str());
    if (verbose) {
        auto bounds = ("(" + std::to_string(lo_error) + ", " + std::to_string(hi_error) + ")").c_str();
        printf("\t↝ search space=%-15s", bounds);
    }
    printf("\n");
}

template<typename K>
void minimize_space_given_time(size_t max_time, float tolerance, std::vector<K> &data,
                               size_t lo_error, size_t hi_error, bool verbose) {
    auto latency = 82.1;
    size_t cache_line = x86_cache_line();
    auto block_size = cache_line / sizeof(int64_t);
    auto starting_error = std::clamp(size_t(block_size * std::pow(2., max_time / latency - 1.)), lo_error, hi_error);
    std::vector<IndexStats> all_stats;

    const size_t starting_i = 2048;
    size_t i = starting_i;
    size_t bin_search_lo = lo_error;
    size_t bin_search_hi = hi_error;
    all_stats.emplace_back(data, starting_error);
    minimize_time_logging(all_stats.back(), verbose, starting_error, starting_error);

    if (all_stats.back().lookup_time < max_time) {
        while (starting_error + (i << 1) < hi_error && all_stats.back().lookup_time < max_time * (1 + tolerance)) {
            i <<= 1;
            all_stats.emplace_back(data, starting_error + i);
            bin_search_lo = starting_error + i / 2;
            bin_search_hi = starting_error + i;
            minimize_time_logging(all_stats.back(), verbose, bin_search_lo, bin_search_hi);
        }
        bin_search_lo = i <= (starting_i << 1) ? starting_error : starting_error + i / 2;
    } else {
        while (starting_error > (i << 1) + lo_error && all_stats.back().lookup_time > max_time * (1 - tolerance)) {
            i <<= 1;
            all_stats.emplace_back(data, starting_error - i);
            bin_search_lo = starting_error - i;
            bin_search_hi = starting_error - i / 2;
            minimize_time_logging(all_stats.back(), verbose, bin_search_lo, bin_search_hi);
        }
        bin_search_hi = i <= (starting_i << 1) ? starting_error : starting_error - i / 2;
    }

    while (bin_search_hi - bin_search_lo > cache_line / 2) {
        i = (bin_search_hi + bin_search_lo) / 2;
        all_stats.emplace_back(data, i);
        if (all_stats.back().lookup_time > max_time)
            bin_search_hi = i;
        else
            bin_search_lo = i + 1;
        minimize_time_logging(all_stats.back(), verbose, bin_search_lo, bin_search_hi);
    }

    auto predicate = [max_time, tolerance](const IndexStats &a) {
        return a.lookup_time <= max_time * (1 + tolerance);
    };
    if (!std::any_of(all_stats.cbegin(), all_stats.cend(), predicate)) {
        printf("It is not possible to satisfy the given constraint. Increase the maximum time.");
        exit(1);
    }

    auto comp = [max_time, tolerance, predicate](const IndexStats &a, const IndexStats &b) {
        return !predicate(b) || (predicate(a) && a.size_in_bytes < b.size_in_bytes);
    };
    auto best = std::min_element(all_stats.cbegin(), all_stats.cend(), comp);

    printf("%s\n", std::string(80, '-').c_str());
    auto time = std::accumulate(all_stats.cbegin(), all_stats.cend(), 0ul,
                                [](size_t result, const IndexStats &s) { return result + s.construction_time; });
    printf("%zu iterations for a total construction time of %.0f s\n", all_stats.size(), time * 1.e-9);
    printf("Set the error to %zu for an index of %zu bytes\n", best->error, best->size_in_bytes);
}

template<typename K>
void minimize_time_given_space(size_t max_space, float tolerance, std::vector<K> &data,
                               size_t lo_error, size_t hi_error, bool verbose) {
    const auto guess_steps_threshold = size_t(2 * std::log2(std::log2(hi_error - lo_error)));
    size_t guess_steps = 0;
    std::vector<IndexStats> all_index_stats;

    double p[2] = {data.size() / 2., 1.000000001};
    alglib::real_1d_array c_starting_point;
    c_starting_point.setcontent(2, p);
    alglib::real_1d_array c_space = c_starting_point;

    do {
        size_t guess = 0;
        size_t mid = (lo_error + hi_error) / 2;

        if (all_index_stats.size() >= 4 && guess_steps < guess_steps_threshold) {
            auto[c_fit, fit_report] = fit_segments_count_model(all_index_stats, c_starting_point);
            if (fit_report.r2 > 0.8) // use the guess of the position only if the fitted model is good
                c_space = c_fit;
        }

        if (guess_steps < guess_steps_threshold) {
            auto constants = sizeof(K) + 2 * sizeof(double);

            guess = size_t(guess_error_space(100, c_space[0], c_space[1], max_space, constants));
            guess = std::clamp(guess, lo_error + 1, hi_error - 1);

            auto bias_weight = guess_steps <= 1 ? 0 : float(guess_steps) / guess_steps_threshold;
            auto biased_guess = mid * bias_weight + guess * (1 - bias_weight);

            mid = size_t(biased_guess);
            guess_steps++;
        }

        IndexStats stats(data, mid);
        all_index_stats.push_back(stats);
        auto kib = stats.size_in_bytes / double(1u << 10u);
        auto query_time = std::to_string(stats.lookup_time) + "±" + std::to_string(stats.lookup_time_std);
        printf("%-19zu %-19.1f %-19.2f %-19s", mid, stats.construction_time * 1.e-9, kib, query_time.c_str());
        if (verbose) {
            auto bounds = ("(" + std::to_string(lo_error) + ", " + std::to_string(hi_error) + ")").c_str();
            printf("\t↝ search space=%-15s \ts(ε)=%.0fε^-%.2f \tε guess=%zu", bounds, c_space[0], c_space[1], guess);
        }
        printf("\n");
        fflush(stdout);

        if (stats.size_in_bytes <= max_space)
            hi_error = mid;
        else
            lo_error = mid + 1;
    } while (lo_error < hi_error
        && std::fabs(all_index_stats.back().size_in_bytes - max_space) > max_space * tolerance);

    printf("%s\n", std::string(80, '-').c_str());
    auto time = std::accumulate(all_index_stats.cbegin(), all_index_stats.cend(), 0ul,
                                [](size_t result, const IndexStats &s) { return result + s.construction_time; });
    printf("Total construction time %.0f s\n", time * 1.e-9);
    printf("Set the error to %zu\n", lo_error);
}