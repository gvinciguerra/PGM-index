/*
 * This example shows how to use pgm::MappedPGMIndex, a container that stores its data on disk and uses the PGM-index to
 * support efficient query operations.
 * Compile with:
 *   g++ mapped.cpp -std=c++17 -I../include -o mapped
 * Run with:
 *   ./mapped
 */

#include <vector>
#include <iterator>
#include "pgm/pgm_index_variants.hpp"

template<typename T, typename Iter>
void print_iterator(Iter first, Iter last) {
    std::cout << "[";
    if (first != last) {
        std::copy(first, std::prev(last), std::ostream_iterator<T>(std::cout, " "));
        std::cout << *std::prev(last);
    }
    std::cout << "]" << std::endl;
}

int main() {
    std::string tmp_filename = std::tmpnam(nullptr);

    {
        // Generate some random data
        std::vector<int> data(1000000);
        std::generate(data.begin(), data.end(), []() { return std::rand() % 100000; });
        std::sort(data.begin(), data.end());

        // Construct the disk-backed container
        pgm::MappedPGMIndex<int, 32> pgm(data.begin(), data.end(), tmp_filename);
        std::cout << "Created indexed file at " << tmp_filename << std::endl
                  << "- elements " << pgm.size() << std::endl
                  << "- file size " << pgm.file_size_in_bytes() << " bytes" << std::endl
                  << "- index size " << pgm.size_in_bytes() << " bytes" << std::endl << std::endl;

        // Here, the container is ready to use
        // When the code block ends, the data and pgm objects gets deallocated
    }

    {
        // Read the existing disk-backed container
        pgm::MappedPGMIndex<int, 32> pgm(tmp_filename);
        std::cout << "File re-opened" << std::endl
                  << "count(15) = " << pgm.count(15) << std::endl
                  << "contains(1) = " << (pgm.contains(1) ? "true" : "false") << std::endl
                  << "lower_bound(42) = " << *pgm.lower_bound(42) << std::endl
                  << "upper_bound(70) = " << *pgm.upper_bound(70) << std::endl
                  << "Range search [50, 60] = ";
        print_iterator<int>(pgm.lower_bound(50), pgm.upper_bound(60));
    }

    std::remove(tmp_filename.c_str());
    return 0;
}