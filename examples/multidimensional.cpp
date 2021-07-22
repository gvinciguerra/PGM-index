/*
 * This example shows how to use pgm::MultidimensionalPGMIndex, a container supporting orthogonal range queries
 * in k dimensions. To run it, your CPU must support the BMI2 instruction set.
 * Compile with:
 *   g++ multidimensional.cpp -std=c++17 -I../include -o multidimensional
 * Run with:
 *   ./multidimensional
 */

#include <vector>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "pgm/pgm_index_variants.hpp"

int main() {
    // Generate random points in a 3D space
    auto rand_coord = []() { return std::rand() % 500; };
    std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> data(1000000);
    std::generate(data.begin(), data.end(), [&] { return std::make_tuple(rand_coord(), rand_coord(), rand_coord()); });

    // Construct the Multidimensional PGM-index
    constexpr auto dimensions = 3; // number of dimensions
    constexpr auto epsilon = 32;   // space-time trade-off parameter
    pgm::MultidimensionalPGMIndex<dimensions, uint64_t, epsilon> pgm_3d(data.begin(), data.end());

    // Query the Multidimensional PGM-index
    std::cout << "range({4,1,2}, {15,20,10}) = ";
    for (auto it = pgm_3d.range({4, 1, 2}, {15, 20, 10}); it != pgm_3d.end(); ++it)
        std::cout << "(" << std::get<0>(*it) << "," << std::get<1>(*it) << "," << std::get<2>(*it) << ") ";
    std::cout << std::endl;

    auto count = std::distance(pgm_3d.range({0, 0, 0}, {5, 10, 15}), pgm_3d.end());
    std::cout << "points in range({0,0,0}, {5,10,15}) = " << count << std::endl;

    auto knn = pgm_3d.knn({4, 1, 2}, 5);
    std::cout << "5 nearest points from {4,1,2} = ";
    for (auto point : knn){
        std::cout << "(" << std::get<0>(point) << "," << std::get<1>(point) << "," << std::get<2>(point) << ") ";
    }
    std::cout << std::endl;
    return 0;
}