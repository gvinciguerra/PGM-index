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

#include <cmath>
#include <limits>
#include <vector>
#include <cassert>
#include <stdexcept>
#include <type_traits>
#include "utils.hpp"

template<typename T>
using LargeSigned = typename std::conditional<std::is_floating_point<T>::value,
                                              long double,
                                              typename std::conditional<(sizeof(T) < 8),
                                                                        int64_t,
                                                                        __int128>::type>::type;

template<typename X, typename Y, typename Floating = double>
class OptimalPiecewiseLinearModel {
private:
    using SX = LargeSigned<X>;
    using SY = LargeSigned<Y>;

    struct Point {
        SX x{};
        SY y{};

        Point() = default;

        Point(SX x, SY y) : x(x), y(y) {};

        Point operator-(Point p) const {
            return Point{SX(x - p.x), SY(y - p.y)};
        }

        bool operator<(Point p) const {
            return y * p.x < x * p.y;
        }

        bool operator>(Point p) const {
            return y * p.x > x * p.y;
        }

        bool operator==(Point p) const {
            return y * p.x == x * p.y;
        }
    };

    const SY error_fwd;
    const SY error_bwd;
    std::vector<Point> lower;
    std::vector<Point> upper;
    size_t lower_start = 0;
    size_t upper_start = 0;
    size_t points_in_hull = 0;
    Point rectangle[4];

    template<typename P>
    SX cross(const P &O, const P &A, const P &B) {
        return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
    }

public:
    explicit OptimalPiecewiseLinearModel(size_t error_fwd, size_t error_bwd)
        : error_fwd(error_fwd), error_bwd(error_bwd) {
        upper.reserve(1u << 16);
        lower.reserve(1u << 16);
    }

    bool add_point(X x, Y y) {
        if (x < rectangle[2].x || x < rectangle[3].x)
            throw std::logic_error("Points must be increasing by x.");

        SX xx = x;
        SY yy = y;
        
        if (points_in_hull == 0) {
            rectangle[0] = {xx, yy + error_fwd};
            rectangle[1] = {xx, yy - error_bwd};
            ++points_in_hull;
            return true;
        }

        if (points_in_hull == 1) {
            rectangle[2] = {xx, yy - error_bwd};
            rectangle[3] = {xx, yy + error_fwd};
            upper.clear();
            upper.push_back(rectangle[0]);
            upper.push_back(rectangle[3]);
            lower.clear();
            lower.push_back(rectangle[1]);
            lower.push_back(rectangle[2]);
            upper_start = lower_start = 0;
            ++points_in_hull;
            return true;
        }

        Point p1(xx, yy + error_fwd);
        Point p2(xx, yy - error_bwd);
        auto slope1 = rectangle[2] - rectangle[0];
        auto slope2 = rectangle[3] - rectangle[1];
        bool outside_line1 = p1 - rectangle[2] < slope1;
        bool outside_line2 = p2 - rectangle[3] > slope2;

        if (outside_line1 || outside_line2) {
            points_in_hull = 0;
            return false;
        }

        if (p1 - rectangle[1] < slope2) {
            // Find extreme slope
            Point min = lower[lower_start] - p1;
            size_t min_i = lower_start;
            for (size_t i = lower_start + 1; i < lower.size(); i++) {
                auto val = (lower[i] - p1);
                if (val > min)
                    break;
                else {
                    min = val;
                    min_i = i;
                }
            }

            rectangle[1] = lower[min_i];
            rectangle[3] = p1;
            lower_start = min_i;

            // Hull update
            size_t end = upper.size();
            for (; end >= upper_start + 2 && cross(upper[end - 2], upper[end - 1], p1) <= 0; --end);
            upper.resize(end);
            upper.push_back(p1);
        }

        if (p2 - rectangle[0] > slope1) {
            // Find extreme slope
            Point max = upper[upper_start] - p2;
            size_t max_i = upper_start;
            for (size_t i = upper_start + 1; i < upper.size(); i++) {
                auto val = (upper[i] - p2);
                if (val < max)
                    break;
                else {
                    max = val;
                    max_i = i;
                }
            }

            rectangle[0] = upper[max_i];
            rectangle[2] = p2;
            upper_start = max_i;

            // Hull update
            size_t end = lower.size();
            for (; end >= lower_start + 2 && cross(lower[end - 2], lower[end - 1], p2) >= 0; --end);
            lower.resize(end);
            lower.push_back(p2);
        }

        ++points_in_hull;
        return true;
    }

    std::pair<double, double> get_intersection() const {
        auto &p0 = rectangle[0];
        auto &p1 = rectangle[1];
        auto &p2 = rectangle[2];
        auto &p3 = rectangle[3];
        auto slope1 = p2 - p0;
        auto slope2 = p3 - p1;

        if (points_in_hull == 1 || slope1 == slope2)
            return std::make_pair(double(p0.x), double(p0.y));

        double a = slope1.x * slope2.y - slope1.y * slope2.x;
        double b = ((p1.x - p0.x) * (p3.y - p1.y) - (p1.y - p0.y) * (p3.x - p1.x)) / a;
        auto i_x = p0.x + b * slope1.x;
        auto i_y = p0.y + b * slope1.y;
        return std::make_pair(i_x, i_y);
    }

    double get_intercept(X key) const {
        auto[i_x, i_y] = get_intersection();
        auto[min_slope, max_slope] = get_slope_range();
        auto slope = 0.5 * (min_slope + max_slope);
        return i_y - (i_x - key) * slope;
    }

    std::pair<double, double> get_slope_range() const {
        if (points_in_hull == 1)
            return {0, 1};
        auto min_slope = double(rectangle[2].y - rectangle[0].y) / (rectangle[2].x - rectangle[0].x);
        auto max_slope = double(rectangle[3].y - rectangle[1].y) / (rectangle[3].x - rectangle[1].x);
        return {min_slope, max_slope};
    }

    void reset() {
        points_in_hull = 0;
    }
};