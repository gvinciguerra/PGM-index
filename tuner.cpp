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

#include <vector>
#include <chrono>
#include <numeric>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "args.hxx"
#include "tuner.hpp"

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

int main(int argc, char **argv) {
    using namespace args;
    ArgumentParser p("Space-time trade-off tuner for the PGM-index.",
                     "This program lets you specify a maximum space and get the PGM-index minimising the query time "
                     "within that space.  Or, it lets you specify a maximum query time and get the PGM-index "
                     "minimising the space.");
    HelpFlag help(p, "help", "Display this help menu", {'h', "help"});
    Group g(p, "Operation modes:", args::Group::Validators::Xor, args::Options::Required);
    ValueFlag<size_t> time(g, "ns", "Specify a time to minimise the space", {'t', "time"});
    ValueFlag<size_t> space(g, "bytes", "Specify a space to minimise the time", {'s', "space"});
    ValueFlag<float> tol(p, "float", "Tolerance between 0-1 on the target resource (default 0.01)", {'o', "tol"}, 0.01);
    Flag verbose(p, "verbose", "Show additional logging info", {'v', "verbose"});
    Group t(p, "File type:", args::Group::Validators::Xor, args::Options::Required);
    Flag bin(t, "binary", "The input file is a binary file containing 32-bit integers", {'b', "binary"});
    Flag csv(t, "csv", "The input file is a csv file containing integers separated by a newline", {'c', "csv"});
    Positional<std::string> file(p, "file", "The file containing the input data", args::Options::Required);
    CompletionFlag completion(p, {"complete"});

    try {
        p.ParseCLI(argc, argv);
    }
    catch (args::Help &) {
        std::cout << p;
        return 0;
    }
    catch (args::Error &e) {
        std::cerr << p;
        return 1;
    }

    std::vector<int64_t> data;
    try {
        data = bin.Get() ? read_data_binary<int32_t, int64_t>(file.Get()) : read_data_csv(file.Get());
    } catch (std::ios_base::failure &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << std::strerror(errno) << std::endl;
        exit(1);
    }

    std::sort(data.begin(), data.end());
    const size_t lo_error = 2 * x86_cache_line() / sizeof(int64_t);
    const size_t hi_error = data.size() / 2;
    bool minimize_space = time;

    printf("Dataset: %zu entries\n", data.size());
    if (minimize_space)
        printf("Max time: %zu±%.0f ns\n", time.Get(), time.Get() * tol.Get());
    else
        printf("Max space: %zu±%.0f KiB\n", space.Get() / (1ul << 10ul), space.Get() * tol.Get() / (1ul << 10ul));

    printf("%s\n", std::string(80, '-').c_str());
    printf("%-19s %-19s %-19s %-19s\n", "Error", "Construction (s)", "Space (KiB)", "Query (ns)");
    printf("%s\n", std::string(80, '-').c_str());

    if (minimize_space)
        minimize_space_given_time(time.Get(), tol.Get(), data, lo_error, hi_error, verbose.Get());
    else
        minimize_time_given_space(space.Get(), tol.Get(), data, lo_error, hi_error, verbose.Get());
}