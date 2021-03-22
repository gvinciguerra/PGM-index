<p align="center">
  <img src="https://pgm.di.unipi.it/images/logo.svg" alt="The PGM-index" style="width: 565px">
</p>

<p align="center">The Piecewise Geometric Model index (PGM-index) is a data structure that enables fast lookup, predecessor, range searches and updates in arrays of billions of items using orders of magnitude less space than traditional indexes while providing the same worst-case query time guarantees.</p>

<p align="center">
    <a href="https://pgm.di.unipi.it/">Website</a>
    | <a href="https://pgm.di.unipi.it/docs">Documentation</a>
    | <a href="http://www.vldb.org/pvldb/vol13/p1162-ferragina.pdf">Paper</a>
    | <a href="https://pgm.di.unipi.it/slides-pgm-index-vldb.pdf">Slides</a>
    | <a href="https://github.com/gvinciguerra/PyGM">Python wrapper</a>
    | <a href="http://acube.di.unipi.it">A³ Lab</a>
</p>

<p align="center">
    <a href="https://github.com/gvinciguerra/PGM-index/actions?query=workflow%3Abuild"><img src="https://img.shields.io/github/workflow/status/gvinciguerra/PGM-index/build" alt="GitHub Workflow Status"></a>
    <a href="https://lgtm.com/projects/g/gvinciguerra/PGM-index/context:cpp"><img alt="Language grade: C/C++" src="https://img.shields.io/lgtm/grade/cpp/github/gvinciguerra/PGM-index?label=code%20quality"/></a>
    <a href="https://github.com/gvinciguerra/PGM-index/blob/master/LICENSE"><img src="https://img.shields.io/github/license/gvinciguerra/PGM-index" alt="License"></a>
    <a href="https://github.com/gvinciguerra/PGM-index/stargazers"><img src="https://img.shields.io/github/stars/gvinciguerra/PGM-index" alt="GitHub stars"></a>
    <a href="https://github.com/gvinciguerra/PGM-index/network/members"><img alt="GitHub forks" src="https://img.shields.io/github/forks/gvinciguerra/PGM-index"></a>
    <a href="https://repl.it/github/gvinciguerra/PGM-index"><img alt="Run on Repl.it" src="https://img.shields.io/badge/run-examples-667881?logo=repl.it&logoColor=white"></a>
</p>

## Quickstart

This is a header-only library. It does not need to be installed. Just clone the repo with

```bash
git clone https://github.com/gvinciguerra/PGM-index.git
cd PGM-index
```

and copy the `include/pgm` directory to your system's or project's include path.
                                                                          
The `examples/simple.cpp` file shows how to index and query a vector of random integers with the PGM-index: 

```cpp
#include <vector>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "pgm/pgm_index.hpp"

int main() {
    // Generate some random data
    std::vector<int> data(1000000);
    std::generate(data.begin(), data.end(), std::rand);
    data.push_back(42);
    std::sort(data.begin(), data.end());

    // Construct the PGM-index
    const int epsilon = 128; // space-time trade-off parameter
    pgm::PGMIndex<int, epsilon> index(data);

    // Query the PGM-index
    auto q = 42;
    auto range = index.search(q);
    auto lo = data.begin() + range.lo;
    auto hi = data.begin() + range.hi;
    std::cout << *std::lower_bound(lo, hi, q);

    return 0;
}
```

[Run and edit this and other examples on Repl.it](https://repl.it/github/gvinciguerra/PGM-index). Or run it locally via:

```bash
g++ examples/simple.cpp -std=c++17 -I./include -o simple
./simple
```

## Classes overview

Other than the `pgm::PGMIndex` class in the example above, this library provides the following classes:

- `pgm::DynamicPGMIndex` supports insertions and deletions.
- `pgm::MultidimensionalPGMIndex` stores points in k dimensions and supports orthogonal range queries. 
- `pgm::MappedPGMIndex` stores data on disk and uses a PGMIndex for fast search operations.
- `pgm::CompressedPGMIndex` compresses the segments to reduce the space usage of the index.
- `pgm::OneLevelPGMIndex` uses a binary search on the segments rather than a recursive structure.
- `pgm::BucketingPGMIndex` uses a top-level lookup table to speed up the search on the segments. 
- `pgm::EliasFanoPGMIndex` uses a top-level succinct structure to speed up the search on the segments.

The full documentation is available [here](https://pgm.di.unipi.it/docs/).

## Compile the tests and the tuner

After cloning the repository, build the project with

```bash
cmake . -DCMAKE_BUILD_TYPE=Release
make -j8
```

The test runner will be placed in `test/`. The [tuner](https://pgm.di.unipi.it/docs/tuner/) executable will be placed in `tuner/`. The [benchmark](https://pgm.di.unipi.it/docs/benchmark/) executable will be placed in `benchmark/`.

## License

This project is licensed under the terms of the Apache License 2.0.

If you use the library please put a link to the [website](https://pgm.di.unipi.it) and cite the following paper:

> Paolo Ferragina and Giorgio Vinciguerra. The PGM-index: a fully-dynamic compressed learned index with provable worst-case bounds. PVLDB, 13(8): 1162-1175, 2020.

```tex
@article{Ferragina:2020pgm,
  Author = {Paolo Ferragina and Giorgio Vinciguerra},
  Title = {The {PGM-index}: a fully-dynamic compressed learned index with provable worst-case bounds},
  Year = {2020},
  Volume = {13},
  Number = {8},
  Pages = {1162--1175},
  Doi = {10.14778/3389133.3389135},
  Url = {https://pgm.di.unipi.it},
  Issn = {2150-8097},
  Journal = {{PVLDB}}}
```
