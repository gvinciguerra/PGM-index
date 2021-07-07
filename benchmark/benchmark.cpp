// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2021 Giorgio Vinciguerra.
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

#include "benchmark.hpp"
#include "args.hxx"
#include "pgm/pgm_index.hpp"
#include "pgm/pgm_index_variants.hpp"

#include <fstream>
#include <functional>
#include <utility>

#define FOR_EACH_EPS(C, K) C<K, 8>, C<K, 16>, C<K, 32>, C<K, 64>, C<K, 128>, C<K, 256>,  C<K, 512>, C<K, 1024>
#define FOR_EACH_EPS_2(C, K, D) C<K, 8, D>, C<K, 16, D>, C<K, 32, D>, C<K, 64, D>, C<K, 128, D>, C<K, 256, D>, \
                                C<K, 512, D>, C<K, 1024, D>

#define FOR_EACH_BPGM(C, K) FOR_EACH_EPS_2(C, K, 1 << 16), FOR_EACH_EPS_2(C, K, 1 << 20), FOR_EACH_EPS_2(C, K, 1 << 24)

#define PGM_CLASSES(K) FOR_EACH_EPS(pgm::PGMIndex, K)
#define BPGM_CLASSES(K) FOR_EACH_BPGM(pgm::BucketingPGMIndex, K)
#define EFPGM_CLASSES(K) FOR_EACH_EPS(pgm::EliasFanoPGMIndex, K)
#define CPGM_CLASSES(K) FOR_EACH_EPS(pgm::CompressedPGMIndex, K)

#define ALL_CLASSES(K) PGM_CLASSES(K), BPGM_CLASSES(K), EFPGM_CLASSES(K), CPGM_CLASSES(K)

template<typename K>
void read_ints_helper(args::PositionalList<std::string> &files,
                      size_t record_size,
                      double lookup_ratio,
                      const std::string &workload) {
    OUT_VERBOSE("Running with " << sizeof(K) << "-byte keys + " << record_size - sizeof(K) << "-byte values")
    for (const auto &file : files.Get()) {
        auto data = to_records(read_data_binary<K>(file, true), record_size);
        auto filename = file.substr(file.find_last_of("/\\") + 1);
        benchmark_all<K, ALL_CLASSES(K)>(filename, data, record_size, lookup_ratio, workload);
    }
}

int main(int argc, char **argv) {
    using namespace args;
    ArgumentParser p("Benchmark for the PGM-index library.");
    p.helpParams.flagindent = 2;
    p.helpParams.helpindent = 25;
    p.helpParams.progindent = 0;
    p.helpParams.descriptionindent = 0;

    CompletionFlag completion(p, {"complete"});
    HelpFlag help(p, "help", "Display this help menu", {'h', "help"});
    Flag verbose(p, "", "Verbose output", {'v', "verbose"});
    ValueFlag<size_t> value_size(p, "bytes", "Size of the values associated to keys", {'V', "values"}, 0);

    Group g1(p, "QUERY WORKLOAD OPTIONS (mutually exclusive):", Group::Validators::AtMostOne);
    ValueFlag<double> ratio(g1, "ratio", "Random workload with the given lookup ratio", {'r', "ratio"}, 0.333);
    ValueFlag<std::string> workload(g1, "file", "Custom workload file. Obeys the format of input files", {'w', "workload"});

    Group g2(p, "INPUT DATA OPTIONS (mutually exclusive):", Group::Validators::Xor, Options::Required);
    ValueFlag<size_t> synthetic(g2, "size", "Generate synthetic data of the given size", {'s', "synthetic"}, 100000000);
    Flag u64(g2, "", "Input files contain unsigned 64-bit ints", {'U', "u64"});
    Flag i64(g2, "", "Input files contain signed 64-bit ints", {'I', "i64"});
    PositionalList<std::string> files(p, "file", "The input files");

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
        std::cerr << e.what() << std::endl;
        std::cerr << p;
        return 1;
    }

    if (ratio.Get() < 0 || ratio.Get() > 1) {
        std::cerr << "Argument to --" << ratio.GetMatcher().GetLongOrAny().str() << " must be between 0.0 and 1.0.";
        return 1;
    }

    if (synthetic && synthetic.Get() < 1000) {
        std::cerr << "Argument to --" << synthetic.GetMatcher().GetLongOrAny().str() << " must be greater than 1000.";
        return 1;
    }

    global_verbose = verbose.Get();
    std::cout << "dataset,class_name,build_ms,bytes,query_ns" << std::endl;

    if (synthetic) {
        auto record_size = value_size.Get() + sizeof(uint64_t);
        auto n = synthetic.Get();
        std::mt19937 generator(std::random_device{}());
        auto gen = [&](auto distribution) {
            std::vector<uint64_t> out(n);
            std::generate(out.begin(), out.end(), [&] { return distribution(generator); });
            std::sort(out.begin(), out.end());
            return to_records(out, record_size);
        };
        std::vector<std::pair<std::string, std::function<std::vector<char>()>>> distributions = {
            {"uniform_dense", std::bind(gen, std::uniform_int_distribution<uint64_t>(0, n * 1000))},
            {"uniform_sparse", std::bind(gen, std::uniform_int_distribution<uint64_t>(0, n * n))},
            {"binomial", std::bind(gen, std::binomial_distribution<uint64_t>(1ull << 50))},
            {"negative_binomial", std::bind(gen, std::negative_binomial_distribution<uint64_t>(1ull << 50, 0.3))},
            {"geometric", std::bind(gen, std::geometric_distribution<uint64_t>(1e-10))},
        };
        OUT_VERBOSE("Generating " << to_metric(n) << " elements (8-byte keys + " << value_size.Get() << "-byte values)")
        OUT_VERBOSE("Total memory for data is " << to_metric(n * record_size, 2, true) << "B")
        for (auto&[name, gen_data] : distributions)
            benchmark_all<uint64_t, ALL_CLASSES(uint64_t)>(name, gen_data(), record_size, ratio.Get(), workload.Get());
    }

    if (i64.Get())
        read_ints_helper<int64_t>(files, value_size.Get() + sizeof(int64_t), ratio.Get(), workload.Get());
    if (u64.Get())
        read_ints_helper<uint64_t>(files, value_size.Get() + sizeof(uint64_t), ratio.Get(), workload.Get());

    return 0;
}