#include <stdio.h>
#include <stdlib.h>
#include "cpgm.h"

int cmp(const void *a, const void *b) {
    uint64_t aa = *(uint64_t *) a;
    uint64_t bb = *(uint64_t *) b;
    if (aa < bb)
        return -1;
    if (aa > bb)
        return 1;
    return 0;
}

int main() {
    // Generate some random data
    size_t n = 1000000;
    uint64_t data[n];
    data[0] = 42;
    for (int i = 1; i < n; ++i)
        data[i] = rand() % 10000000;
    qsort(data, n, sizeof(data[0]), cmp);

    // Construct the PGM-index
    const size_t epsilon = 32; // space-time trade-off parameter
    pgm_index_uint64_t *pgm = pgm_index_uint64_create(data, n, epsilon);
    printf("PGM-index takes %zu bytes\n", pgm_index_uint64_size_in_bytes(pgm));

    // Query the PGM-index
    uint64_t q = 42;
    approx_pos_t range = pgm_index_uint64_search(pgm, q);
    uint64_t *ptr = (uint64_t *) bsearch(&q, data + range.lo, range.hi - range.lo, sizeof(data[0]), cmp);
    if (ptr)
        printf("Found %llu at position %zu\n", q, ptr - data);
    else
        printf("Not found\n");

    pgm_index_uint64_destroy(pgm);

    return 0;
}