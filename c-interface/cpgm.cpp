// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2020 Giorgio Vinciguerra.
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

#include "cpgm.h"
#include "pgm/pgm_index.hpp"
#include "pgm/pgm_index_dynamic.hpp"

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <type_traits>
#include <vector>

#define EPSILON_RECURSIVE 4

template<typename K>
class PGMWrapper : public pgm::PGMIndex<K, 1, EPSILON_RECURSIVE> {
    size_t epsilon;

public:

    PGMWrapper(const K *a, size_t n, size_t epsilon) : epsilon(epsilon) {
        this->n = n;
        this->first_key = n ? *a : 0;
        this->build(a, a + n, epsilon, EPSILON_RECURSIVE, this->segments, this->levels_offsets);
    }

    approx_pos_t search(const K &key) const {
        auto k = std::max(this->first_key, key);
        auto it = this->segment_for_key(k);
        auto pos = std::min<size_t>((*it)(k), std::next(it)->intercept);
        auto lo = PGM_SUB_EPS(pos, epsilon);
        auto hi = PGM_ADD_EPS(pos, epsilon, this->n);
        return {pos, lo, hi};
    }
};

#define PGM_INDEX_DEFINE(type)                                                                                         \
    struct pgm_index_##type##_ : public PGMWrapper<PGM_T(type)> {                                                      \
        pgm_index_##type##_(const PGM_T(type) * a, size_t n, size_t epsilon)                                           \
            : PGMWrapper<PGM_T(type)>(a, n, epsilon) {}                                                                \
    };                                                                                                                 \
                                                                                                                       \
    PGM_PTR(pgm_index, type)                                                                                           \
    pgm_index_##type##_create(const PGM_T(type) * a, size_t n, size_t epsilon) {                                       \
        try {                                                                                                          \
            return new pgm_index_##type##_(a, n, epsilon);                                                             \
        } catch (const std::invalid_argument &) {                                                                      \
            return nullptr;                                                                                            \
        }                                                                                                              \
    }                                                                                                                  \
                                                                                                                       \
    void pgm_index_##type##_destroy(PGM_PTR(pgm_index, type) pgm) { delete pgm; }                                      \
                                                                                                                       \
    approx_pos_t pgm_index_##type##_search(PGM_PTR(pgm_index, type) pgm, PGM_T(type) q) { return pgm->search(q); }     \
                                                                                                                       \
    size_t pgm_index_##type##_size_in_bytes(PGM_PTR(pgm_index, type) pgm) { return pgm->size_in_bytes(); }

PGM_INDEX_DEFINE(int32)
PGM_INDEX_DEFINE(int64)
PGM_INDEX_DEFINE(uint32)
PGM_INDEX_DEFINE(uint64)

#define DYNAMIC_PGM_INDEX_DEFINE(type)                                                                                 \
    struct dynamic_pgm_index_##type##_ : public pgm::DynamicPGMIndex<PGM_T(type), PGM_T(type)> {                       \
        using pgm::DynamicPGMIndex<PGM_T(type), PGM_T(type)>::DynamicPGMIndex;                                         \
    };                                                                                                                 \
                                                                                                                       \
    PGM_PTR(dynamic_pgm_index, type) dynamic_pgm_index_##type##_create(const pair_##type##_t *a, size_t n) {           \
        try {                                                                                                          \
            return new dynamic_pgm_index_##type##_(a, a + n);                                                          \
        } catch (const std::invalid_argument &) {                                                                      \
            return nullptr;                                                                                            \
        }                                                                                                              \
    }                                                                                                                  \
                                                                                                                       \
    PGM_PTR(dynamic_pgm_index, type) dynamic_pgm_index_##type##_create_empty() {                                       \
        return new dynamic_pgm_index_##type##_();                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    void dynamic_pgm_index_##type##_destroy(PGM_PTR(dynamic_pgm_index, type) pgm) { delete pgm; }                      \
                                                                                                                       \
    size_t dynamic_pgm_index_##type##_size(PGM_PTR(dynamic_pgm_index, type) pgm) { return pgm->size(); }               \
                                                                                                                       \
    size_t dynamic_pgm_index_##type##_size_in_bytes(PGM_PTR(dynamic_pgm_index, type) pgm) {                            \
        return pgm->size_in_bytes();                                                                                   \
    }                                                                                                                  \
                                                                                                                       \
    size_t dynamic_pgm_index_##type##_index_size_in_bytes(PGM_PTR(dynamic_pgm_index, type) pgm) {                      \
        return pgm->index_size_in_bytes();                                                                             \
    }                                                                                                                  \
                                                                                                                       \
    void dynamic_pgm_index_##type##_insert_or_assign(PGM_PTR(dynamic_pgm_index, type) pgm, PGM_T(type) key,            \
                                                     PGM_T(type) value) {                                              \
        pgm->insert_or_assign(key, value);                                                                             \
    }                                                                                                                  \
                                                                                                                       \
    void dynamic_pgm_index_##type##_erase(PGM_PTR(dynamic_pgm_index, type) pgm, PGM_T(type) key) { pgm->erase(key); }  \
                                                                                                                       \
    bool dynamic_pgm_index_##type##_find(PGM_PTR(dynamic_pgm_index, type) pgm, PGM_T(type) key, PGM_T(type) * value) { \
        auto it = pgm->find(key);                                                                                      \
        if (it != pgm->end()) {                                                                                        \
            *value = it->second;                                                                                       \
            return true;                                                                                               \
        }                                                                                                              \
        return false;                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    dynamic_pgm_index_##type##_iterator_t dynamic_pgm_index_##type##_begin(PGM_PTR(dynamic_pgm_index, type) pgm) {     \
        return new dynamic_pgm_index_##type##_::iterator(pgm->begin());                                                \
    }                                                                                                                  \
                                                                                                                       \
    void *dynamic_pgm_index_##type##_lower_bound(PGM_PTR(dynamic_pgm_index, type) pgm, PGM_T(type) q) {                \
        return new dynamic_pgm_index_##type##_::iterator(pgm->lower_bound(q));                                         \
    }                                                                                                                  \
                                                                                                                       \
    bool dynamic_pgm_index_##type##_iterator_next(PGM_PTR(dynamic_pgm_index, type) pgm,                                \
                                                  dynamic_pgm_index_##type##_iterator_t it, PGM_T(type) * key,         \
                                                  PGM_T(type) * value) {                                               \
        auto iter = static_cast<dynamic_pgm_index_##type##_::iterator *>(it);                                          \
        if (pgm->end() == *iter)                                                                                       \
            return false;                                                                                              \
        *key = (*iter)->first;                                                                                         \
        *value = (*iter)->second;                                                                                      \
        ++(*iter);                                                                                                     \
        return true;                                                                                                   \
    }                                                                                                                  \
                                                                                                                       \
    void dynamic_pgm_index_##type##_iterator_destroy(dynamic_pgm_index_##type##_iterator_t it) {                       \
        delete static_cast<dynamic_pgm_index_##type##_::iterator *>(it);                                               \
    }

DYNAMIC_PGM_INDEX_DEFINE(int32)
DYNAMIC_PGM_INDEX_DEFINE(int64)
DYNAMIC_PGM_INDEX_DEFINE(uint32)
DYNAMIC_PGM_INDEX_DEFINE(uint64)