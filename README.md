<p align="center">
  <img src="https://gvinciguerra.github.io/PGM-index/images/logo.svg" alt="The PGM-index" style="width: 565px">
</p>

<p align="center">The Piecewise Geometric Model index (PGM-index) is a data structure that enables fast lookup, predecessor, range searches and updates in arrays of billions of items using orders of magnitude less space than traditional indexes while providing the same worst-case query time guarantees.</p>

<p align="center">
    <a href="https://pgm.di.unipi.it/">Website</a>
    | <a href="https://pgm.di.unipi.it/docs">Documentation</a>
    | <a href="http://www.vldb.org/pvldb/vol13/p1162-ferragina.pdf">Paper</a>
    | <a href="https://github.com/gvinciguerra/PyGM">Python wrapper</a>
    | <a href="http://acube.di.unipi.it">A³ Lab</a>
</p>

<p align="center">
    <a href="https://travis-ci.org/gvinciguerra/PGM-index"><img src="https://img.shields.io/travis/gvinciguerra/PGM-index" alt="Travis (.org)"></a>
    <a href="https://github.com/gvinciguerra/PGM-index/blob/master/LICENSE"><img src="https://img.shields.io/github/license/gvinciguerra/PGM-index" alt="License"></a>
    <a href="https://github.com/gvinciguerra/PGM-index/stargazers"><img src="https://img.shields.io/github/stars/gvinciguerra/PGM-index" alt="GitHub stars"></a>
    <a href="https://github.com/gvinciguerra/PGM-index/network/members"><img alt="GitHub forks" src="https://img.shields.io/github/forks/gvinciguerra/PGM-index"></a>
</p>

## Building the code

To download and build the library use the following commands:

```bash
git clone https://github.com/gvinciguerra/PGM-index.git
cd PGM-index
cmake . -DCMAKE_BUILD_TYPE=Release
make -j8
```

Now you can run the unit tests via:

```
./test/tests
```

## Minimal example

```cpp
#include <vector>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "pgm_index.hpp"

int main(int argc, char **argv) {
    // Generate some random data
    std::vector<int> dataset(1000000);
    std::generate(dataset.begin(), dataset.end(), std::rand);
    dataset.push_back(42);
    std::sort(dataset.begin(), dataset.end());

    // Construct the PGM-index
    const int epsilon = 128; // space-time trade-off parameter
    PGMIndex<int, epsilon> index(dataset);

    // Query the PGM-index
    auto q = 42;
    auto approx_range = index.find_approximate_position(q);
    auto lo = dataset.begin() + approx_range.lo;
    auto hi = dataset.begin() + approx_range.hi;
    std::cout << *std::lower_bound(lo, hi, q);

    return 0;
}
```

## License

This project is licensed under the terms of the GNU General Public License v3.0.

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
