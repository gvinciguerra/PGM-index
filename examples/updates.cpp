/*
 * This example shows how to use pgm::DynamicPGMIndex, a std::map-like container supporting inserts and deletes.
 * Compile with:
 *   g++ updates.cpp -std=c++17 -I../include -o updates
 * Run with:
 *   ./updates
 */

#include <vector>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "pgm/pgm_index_dynamic.hpp"

int main() {
    // Generate some random key-value pairs to bulk-load the Dynamic PGM-index
    std::vector<std::pair<uint32_t, uint32_t>> data(1000000);
    std::generate(data.begin(), data.end(), [] { return std::make_pair(std::rand(), std::rand()); });
    std::sort(data.begin(), data.end());

    // Construct and bulk-load the Dynamic PGM-index
    pgm::DynamicPGMIndex<uint32_t, uint32_t> dynamic_pgm(data.begin(), data.end());

    // Insert some data
    dynamic_pgm.insert_or_assign(2, 4);
    dynamic_pgm.insert_or_assign(4, 8);
    dynamic_pgm.insert_or_assign(8, 16);

    // Delete data
    dynamic_pgm.erase(4);

    // Query the container
    std::cout << "Container size (data + index) = " << dynamic_pgm.size_in_bytes() << " bytes" << std::endl;
    std::cout << "find(4) = " << (dynamic_pgm.find(4) == dynamic_pgm.end() ? "not found" : "found") << std::endl;
    std::cout << "find(8)->second = " << dynamic_pgm.find(8)->second << std::endl;

    std::cout << "Range search [1, 10000) = ";
    auto result = dynamic_pgm.range(1, 10000);
    for (auto[k, v] : result)
        std::cout << "(" << k << "," << v << "), ";

    return 0;
}