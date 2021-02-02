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

#pragma once

#ifdef __cplusplus
#include <cstdint>
#include <cstddef>
extern "C" {
#else
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#endif

typedef struct { size_t pos; size_t lo; size_t hi; } approx_pos_t;

#define PGM_T(x) x##_t
#define PGM_PTR(name, x) name##_##x##_t *

#define PGM_INDEX_DECLARE(type)                                                                                        \
    typedef struct pgm_index_##type##_ pgm_index_##type##_t;                                                           \
    PGM_PTR(pgm_index, type) pgm_index_##type##_create(const PGM_T(type) * a, size_t n, size_t epsilon);               \
    void pgm_index_##type##_destroy(PGM_PTR(pgm_index, type) pgm);                                                     \
    approx_pos_t pgm_index_##type##_search(PGM_PTR(pgm_index, type) pgm, PGM_T(type) q);                               \
    size_t pgm_index_##type##_size_in_bytes(PGM_PTR(pgm_index, type) pgm);

PGM_INDEX_DECLARE(int32)
PGM_INDEX_DECLARE(int64)
PGM_INDEX_DECLARE(uint32)
PGM_INDEX_DECLARE(uint64)

#define DYNAMIC_PGM_INDEX_DECLARE(type)                                                                                \
    typedef struct {                                                                                                   \
        PGM_T(type) first;                                                                                             \
        PGM_T(type) second;                                                                                            \
    } pair_##type##_t;                                                                                                 \
    typedef struct dynamic_pgm_index_##type##_ dynamic_pgm_index_##type##_t;                                           \
    typedef void *dynamic_pgm_index_##type##_iterator_t;                                                               \
    PGM_PTR(dynamic_pgm_index, type) dynamic_pgm_index_##type##_create(const pair_##type##_t *a, size_t n);            \
    PGM_PTR(dynamic_pgm_index, type) dynamic_pgm_index_##type##_create_empty();                                        \
    void dynamic_pgm_index_##type##_destroy(PGM_PTR(dynamic_pgm_index, type) pgm);                                     \
    size_t dynamic_pgm_index_##type##_size(PGM_PTR(dynamic_pgm_index, type) pgm);                                      \
    size_t dynamic_pgm_index_##type##_size_in_bytes(PGM_PTR(dynamic_pgm_index, type) pgm);                             \
    size_t dynamic_pgm_index_##type##_index_size_in_bytes(PGM_PTR(dynamic_pgm_index, type) pgm);                       \
    void dynamic_pgm_index_##type##_insert_or_assign(PGM_PTR(dynamic_pgm_index, type) pgm, PGM_T(type) key,            \
                                                     PGM_T(type) value);                                               \
    void dynamic_pgm_index_##type##_erase(PGM_PTR(dynamic_pgm_index, type) pgm, PGM_T(type) key);                      \
    bool dynamic_pgm_index_##type##_find(PGM_PTR(dynamic_pgm_index, type) pgm, PGM_T(type) key, PGM_T(type) * value);  \
    dynamic_pgm_index_##type##_iterator_t dynamic_pgm_index_##type##_begin(PGM_PTR(dynamic_pgm_index, type) pgm);      \
    dynamic_pgm_index_##type##_iterator_t dynamic_pgm_index_##type##_lower_bound(PGM_PTR(dynamic_pgm_index, type) pgm, \
                                                                                 PGM_T(type) q);                       \
    bool dynamic_pgm_index_##type##_iterator_next(PGM_PTR(dynamic_pgm_index, type) pgm,                                \
                                                  dynamic_pgm_index_##type##_iterator_t it, PGM_T(type) * key,         \
                                                  PGM_T(type) * value);                                                \
    void dynamic_pgm_index_##type##_iterator_destroy(dynamic_pgm_index_##type##_iterator_t it);

DYNAMIC_PGM_INDEX_DECLARE(int32)
DYNAMIC_PGM_INDEX_DECLARE(int64)
DYNAMIC_PGM_INDEX_DECLARE(uint32)
DYNAMIC_PGM_INDEX_DECLARE(uint64)

#ifdef __cplusplus
}
#endif