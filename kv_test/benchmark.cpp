/*
 * Benchmark from existing mapped pgm file
 *
 * This source is derived from examples/mapped.cpp
 */

#include <chrono>
#include <iterator>
#include <math.h>
#include <vector>
#include "pgm/pgm_index_variants.hpp"

#include "flags.h"

#ifndef __has_include
  static_assert(false, "__has_include not supported");
#else
#  if __cplusplus >= 201703L && __has_include(<filesystem>)
#    include <filesystem>
     namespace fs = std::filesystem;
#  elif __has_include(<experimental/filesystem>)
#    include <experimental/filesystem>
     namespace fs = std::experimental::filesystem;
#  elif __has_include(<boost/filesystem.hpp>)
#    include <boost/filesystem.hpp>
     namespace fs = boost::filesystem;
#  endif
#endif

double report_t(size_t t_idx, size_t &count_milestone, size_t &last_count_milestone, long long &last_elapsed, std::chrono::time_point<std::chrono::high_resolution_clock> start_t) {
    const double freq_mul = 1.1;
    auto curr_time = std::chrono::high_resolution_clock::now();
    long long time_elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(curr_time - start_t).count();
    std::cout
        << "t = "<< time_elapsed << " ns: "
        << t_idx + 1 << " counts, tot "
        << (time_elapsed) / (t_idx + 1.0) << "/op, seg "
        << (time_elapsed - last_elapsed) / (t_idx + 1.0 - last_count_milestone) << "/op"
        << std::endl;
    last_elapsed = time_elapsed;
    last_count_milestone = count_milestone;
    count_milestone = ceil(((double) count_milestone) * freq_mul);  // next milestone to print
    return time_elapsed;
}

int main(int argc, char* argv[]) {
    auto flags = parse_flags(argc, argv);

    // extract paths
    std::string key_path = get_required(flags, "key_path");  // query keys and answers
    std::string target_db_path = get_required(flags, "target_db_path");  // path to save db file
    std::string out_path = get_required(flags, "out_path");  // output log
    fs::path target_path(target_db_path);

    // load keyset
    std::vector<uint64_t> queries;
    std::vector<uint64_t> expected_ans;
    {
        std::ifstream query_words_in(key_path);
        std::string line;
        while (std::getline(query_words_in, line)) {
            std::istringstream input;
            input.str(line);

            std::string key;
            std::string exp;
            input >> key;
            input >> exp;

            queries.push_back(std::stoull(key));
            expected_ans.push_back(std::stoull(exp));
        }   
    }

    // variables for milestone
    size_t last_count_milestone = 0;
    size_t count_milestone = 1;
    long long last_elapsed = 0;
    std::vector<double> timestamps;

    // start timer
    auto start_t = std::chrono::high_resolution_clock::now();

    // load mapped pgm
    pgm::MappedPGMIndex<uint64_t, 32, 32> pgm(target_db_path);

    // issue queries and check answers
    size_t err;
    for (size_t t_idx = 0; t_idx < queries.size(); t_idx++) {
        // query and answer
        uint64_t key = queries[t_idx];
        uint64_t answer = expected_ans[t_idx];  

        // search
        auto pgm_idx = pgm.lower_bound(key);
        auto pgm_rank = std::distance(pgm.begin(), pgm_idx);

        // check with answer
        if (pgm_rank != answer) {
            printf("ERROR: incorrect rank: %d, expected: %d\n", pgm_rank, answer);
        }

        if (t_idx + 1 == count_milestone || t_idx + 1 == queries.size()) {
            timestamps.push_back(report_t(t_idx, count_milestone, last_count_milestone, last_elapsed, start_t));    
        }
    }

    // write result to file
    {
        std::cout << "Writing timestamps to file " << out_path << std::endl;
        std::ofstream file_out;
        file_out.open(out_path, std::ios_base::app);
        for (const auto& timestamp : timestamps) {
            file_out << (long long) timestamp << ",";
        }
        file_out << std::endl;
        file_out.close();   
    }
    return 0;
}