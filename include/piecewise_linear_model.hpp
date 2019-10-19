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

#include <list>
#include <cmath>
#include <limits>
#include <vector>
#include <cassert>
#include <stdexcept>
#include <type_traits>
#include "utils.hpp"

template<typename T>
using LargeSigned = typename std::conditional<std::is_integral<T>::value, intmax_t, long double>::type;

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
    };

    const SY error_fwd;
    const SY error_bwd;
    std::list<Point> lower;
    std::list<Point> upper;
    size_t points_in_hull{};
    Point rectangle[4];

    template<typename P>
    SX cross(const P &O, const P &A, const P &B) {
        return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
    }

    auto minmax(std::list<Point> &hull, Point &s) {
        return std::minmax_element(hull.cbegin(), hull.cend(), [&s](auto &u, auto &v) {
            return (u - s) < (v - s);
        });
    }

public:
    explicit OptimalPiecewiseLinearModel(size_t error_fwd, size_t error_bwd)
        : error_fwd(error_fwd), error_bwd(error_bwd) {}

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
            upper = {rectangle[0], rectangle[3]};
            lower = {rectangle[1], rectangle[2]};
            ++points_in_hull;
            return true;
        }

        auto slope1 = rectangle[2] - rectangle[0];
        auto slope2 = rectangle[3] - rectangle[1];

        bool outside_line1 = slope1.y * (xx - rectangle[2].x) > slope1.x * (yy - rectangle[2].y + error_fwd);
        bool outside_line2 = slope2.y * (xx - rectangle[3].x) < slope2.x * (yy - rectangle[3].y - error_bwd);
        if (outside_line1 || outside_line2) {
            points_in_hull = 0;
            return false;
        }

        {
            Point s(xx, yy + error_fwd);
            auto m = minmax(lower, s).first;
            auto new_slope = s - *m;
            if (new_slope < slope2) {
                rectangle[1] = *m;
                rectangle[3] = s;
                lower.erase(lower.begin(), m);
            }

            // Hull update
            while (upper.size() >= 2 && cross(*std::prev(upper.end(), 2), *std::prev(upper.end()), s) <= 0)
                upper.pop_back();
            upper.push_back(s);
        }

        {
            Point s(xx, yy - error_bwd);
            auto m = minmax(upper, s).second;
            auto new_slope = s - *m;
            if (slope1 < new_slope) {
                rectangle[0] = *m;
                rectangle[2] = s;
                upper.erase(upper.begin(), m);
            }

            // Hull update
            while (lower.size() >= 2 && cross(*std::prev(lower.end(), 2), *std::prev(lower.end()), s) >= 0)
                lower.pop_back();
            lower.push_back(s);
        }

        ++points_in_hull;
        return true;
    }

    std::pair<double, double> get_intersection() const {
        auto &p0 = rectangle[0];
        auto &p1 = rectangle[1];
        auto &p2 = rectangle[2];
        auto &p3 = rectangle[3];
        auto slope1 = rectangle[2] - rectangle[0];
        auto slope2 = rectangle[3] - rectangle[1];

        if (points_in_hull == 1 || slope1.x * slope2.y == slope1.y * slope2.x)
            return std::make_pair(double(p0.x), double(p0.y));

        double a = slope1.x * slope2.y - slope1.y * slope2.x;
        double b = ((p1.x - p0.x) * (p3.y - p1.y) - (p1.y - p0.y) * (p3.x - p1.x)) / a;
        auto i_x = rectangle[0].x + b * (rectangle[2].x - rectangle[0].x);
        auto i_y = rectangle[0].y + b * (rectangle[2].y - rectangle[0].y);
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