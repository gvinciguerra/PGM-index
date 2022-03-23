/*
 * Build mapped pgm file
 *
 * This source is derived from examples/mapped.cpp

export dataset_name="wiki_ts_200M_uint64"

export ROOT=$(pwd)/..
mkdir -p ${ROOT}/temp/pgm
mkdir -p ${ROOT}/storage/pgm
mkdir -p ${ROOT}/out

./kv_test/kv_build \
    --data_path=${ROOT}/data/${dataset_name} \
    --build_db_path=${ROOT}/temp/pgm/${dataset_name} \
    --target_db_path=${ROOT}/storage/pgm/${dataset_name}

./kv_test/kv_benchmark \
    --key_path=${ROOT}/keyset/${dataset_name}_ks \
    --target_db_path=${ROOT}/storage/pgm/${dataset_name} \
    --out_path=out/${dataset_name}_out.txt
 */

#include <vector>
#include <iterator>
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


/* Reuse code from benchmark */

bool global_verbose = true;

#define IF_VERBOSE(X) if (global_verbose) { X; }
#define OUT_VERBOSE(X) if (global_verbose) { std::cout << "# " << X << std::endl; }

template<typename T>
std::string to_metric(const T &x, int digits = 2, bool space = false) {
    static const char *prefix[] = {"y", "z", "a", "f", "p", "n", "µ", "m", "", "k", "M", "G", "T", "P", "E", "Z", "Y"};
    double value = x;
    if (value < 0.)
        value = -value;

    auto log_x = (int) std::log10(value);
    if (log_x > 0)
        log_x = (log_x / 3) * 3;
    else
        log_x = (-log_x + 3) / 3 * (-3);

    value *= std::pow(10, -log_x);
    if (value >= 1000.)
        value /= 1000.0, log_x += 3;

    if (std::fmod(value, 1.0) == 0.0)
        digits = 0;

    constexpr auto prefix_start = -24;
    const char *fmt = space ? (x < 0. ? "-%.*f %s" : "%.*f %s") : (x < 0. ? "-%.*f%s" : "%.*f%s");
    const char *chosen_prefix = prefix[(log_x - prefix_start) / 3];
    std::vector<char> buffer(100);
    int size = std::snprintf(&buffer[0], buffer.size(), fmt, digits, value, chosen_prefix);
    return std::string(&buffer[0], size);
}

template<typename T>
std::vector<T> read_data_binary(const std::string &filename, bool check_sorted) {
    OUT_VERBOSE("Reading " << filename)
    std::vector<T> data;
    try {
        std::fstream in(filename, std::ios::in | std::ios::binary);
        in.exceptions(std::ios::failbit | std::ios::badbit);
        struct stat fs;
        stat(filename.c_str(), &fs);
        size_t file_size = fs.st_size;
        data.resize(file_size / sizeof(T));
        in.read((char *) data.data(), file_size);
    } catch (std::ios_base::failure &e) {
        std::cerr << "Could not read the file. " << e.what() << std::endl;
        exit(1);
    }

    if (size_t(data[0]) == data.size() - 1)
        data.erase(data.begin()); // in some formats, the first element is the size of the dataset, ignore it

    if (check_sorted && !std::is_sorted(data.begin(), data.end())) {
        std::cerr << "Input data must be sorted." << std::endl;
        std::cerr << "Read: [";
        std::copy(data.begin(), std::min(data.end(), data.begin() + 10), std::ostream_iterator<T>(std::cerr, ", "));
        std::cout << "...]." << std::endl;
        exit(1);
    }

    IF_VERBOSE(std::cout << "# Read " << to_metric(data.size()) << " elements: [")
    IF_VERBOSE(std::copy(data.begin(),
                         std::min(data.end() - 1, data.begin() + 5),
                         std::ostream_iterator<T>(std::cout, ", ")))
    IF_VERBOSE(std::cout << "..., " << *(data.end() - 1) << "]" << std::endl)

    return data;
}


/* Build! */

int main(int argc, char* argv[]) {
    auto flags = parse_flags(argc, argv);

    // extract paths
    std::string data_path = get_required(flags, "data_path");  // sosd blob
    std::string build_db_path = get_required(flags, "build_db_path");  // path to build (e.g. locally)
    std::string target_db_path = get_required(flags, "target_db_path");  // path to save db file
    fs::path build_path(build_db_path);
    fs::path target_path(target_db_path);

    // Load sosd blob
    std::vector<uint64_t> data = read_data_binary<uint64_t>(data_path, true);

    // Construct the disk-backed container
    pgm::MappedPGMIndex<uint64_t, 32, 32> pgm(data.begin(), data.end(), build_db_path);
    std::cout << "Created indexed file at " << build_db_path << std::endl
              << "- elements " << pgm.size() << std::endl
              << "- file size " << pgm.file_size_in_bytes() << " bytes" << std::endl
              << "- index size " << pgm.size_in_bytes() << " bytes" << std::endl << std::endl;

    // Copy to target directory (slow random write)
    if (build_path != target_path) {
        std::cout << "Copying DB directory from " << build_path << " to " << target_path.string() << std::endl;      
        fs::copy(build_path, target_path, fs::copy_options::overwrite_existing | fs::copy_options::recursive);

        std::cout << "Removing tempporary DB directory " << build_path << std::endl;      
        std::remove(build_db_path.c_str());
    }
    return 0;
}