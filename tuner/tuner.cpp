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

#include "tuner.hpp"
#include "args.hxx"
#include <cstdint>
#include <iostream>

template<typename K>
void run_tuner(args::ValueFlag<size_t> &time,
               args::ValueFlag<size_t> &space,
               args::ValueFlag<double> &tol,
               args::ValueFlag<float> &ratio,
               args::Positional<std::string> &file);

int main(int argc, char **argv) {
    using namespace args;
    ArgumentParser p("Space-time trade-off tuner for the PGM-index. \n\nThis program lets you specify a maximum space "
                     "and get the PGM-index minimising the query time within that space. Or, it lets you specify a "
                     "maximum query time and get the PGM-index minimising the space.");
    p.helpParams.flagindent = 2;
    p.helpParams.helpindent = 25;
    p.helpParams.progindent = 0;
    p.helpParams.descriptionindent = 0;

    HelpFlag help(p, "help", "Display this help menu", {'h', "help"});
    ValueFlag<double> tol(p, "float", "Tolerance between 0 and 1 on the constraint (default 0.01)", {'o', "tol"}, 0.01);
    ValueFlag<float> ratio(p, "ratio", "Ratio of lookups in the query workload (default 0.33)", {'r', "ratio"}, 0.333);
    Flag verbose(p, "verbose", "Show additional logging info", {'v', "verbose"});

    Group g(p, "OPERATION MODES:", args::Group::Validators::Xor, args::Options::Required);
    ValueFlag<size_t> time(g, "ns", "Specify a time to minimise the space", {'t', "time"});
    ValueFlag<size_t> space(g, "bytes", "Specify a space to minimise the time", {'s', "space"});

    Group t(p, "INPUT DATA OPTIONS:", args::Group::Validators::Xor, args::Options::Required);
    Flag u64(t, "", "Input file contains unsigned 64-bit ints", {'U', "u64"});
    Flag i64(t, "", "Input file contains signed 64-bit ints", {'I', "i64"});
    Flag u32(t, "", "Input file contains unsigned 32-bit ints", {'u', "u32"});
    Flag i32(t, "", "Input file contains signed 32-bit ints", {'i', "i32"});
    Positional<std::string> file(p, "file", "The file containing the input data", args::Options::Required);
    CompletionFlag completion(p, {"complete"});

    try {
        p.ParseCLI(argc, argv);
    }
    catch (args::Completion &e) {
        std::cout << e.what();
        return 0;
    }
    catch (args::Help &) {
        std::cout << p;
        return 0;
    }
    catch (args::Error &e) {
        std::cerr << p;
        return 1;
    }

    global_verbose = verbose.Get();

    if (i64.Get())
        run_tuner<int64_t>(time, space, tol, ratio, file);
    if (u64.Get())
        run_tuner<uint64_t>(time, space, tol, ratio, file);
    if (i32.Get())
        run_tuner<int32_t>(time, space, tol, ratio, file);
    if (u32.Get())
        run_tuner<uint32_t>(time, space, tol, ratio, file);
}

template<typename K>
void run_tuner(args::ValueFlag<size_t> &time,
               args::ValueFlag<size_t> &space,
               args::ValueFlag<double> &tol,
               args::ValueFlag<float> &ratio,
               args::Positional<std::string> &file) {
    std::vector<K> data = read_data_binary<K>(file.Get(), true);

    std::sort(data.begin(), data.end());
    auto lo_eps = 2 * cache_line_size() / sizeof(int64_t);
    auto hi_eps = data.size() / 2;
    auto minimize_space = time.Matched();

    std::printf("Dataset: %zu entries\n", data.size());
    if (minimize_space)
        std::printf("Max time: %zu±%.0f ns\n", time.Get(), time.Get() * tol.Get());
    else
        std::printf("Max space: %zu±%.0f KiB\n", space.Get() / (1ul << 10ul), space.Get() * tol.Get() / (1ul << 10ul));

    std::printf("%s\n", std::string(80, '-').c_str());
    std::printf("%-19s %-19s %-19s %-19s\n", "Epsilon", "Construction (s)", "Space (KiB)", "Query (ns)");
    std::printf("%s\n", std::string(80, '-').c_str());

    auto queries = generate_queries(data.begin(), data.end(), ratio.Get(), 1000000);
    if (minimize_space)
        minimize_space_given_time(time.Get(), tol.Get(), data, queries, lo_eps, hi_eps, global_verbose);
    else
        minimize_time_given_space(space.Get(), tol.Get(), data, queries, lo_eps, hi_eps, global_verbose);
}