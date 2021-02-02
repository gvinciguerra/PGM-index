#include <stdio.h>
#include <stdlib.h>
#include "cpgm.h"

int cmp(const void *a, const void *b) {
    const pair_int32_t *aa = (pair_int32_t *) a;
    const pair_int32_t *bb = (pair_int32_t *) b;
    if (aa->first < bb->first)
        return -1;
    if (aa->first > bb->first)
        return 1;
    return 0;
}

int main() {
    // Generate some random key-value pairs to bulk-load the Dynamic PGM-index
    size_t n = 1000000;
    pair_int32_t data[n];
    for (int i = 1; i < n; ++i) {
        data[i].first = rand();
        data[i].second = rand();
    }
    qsort(data, n, sizeof(data[0]), cmp);

    // Construct and bulk-load the Dynamic PGM-index
    dynamic_pgm_index_int32_t *pgm = dynamic_pgm_index_int32_create(data, n);

    // Insert some data
    dynamic_pgm_index_int32_insert_or_assign(pgm, 2, 4);
    dynamic_pgm_index_int32_insert_or_assign(pgm, 4, 8);
    dynamic_pgm_index_int32_insert_or_assign(pgm, 8, 16);

    // Delete data
    dynamic_pgm_index_int32_erase(pgm, 4);

    // Query the container
    int32_t key;
    int32_t value;

    if (dynamic_pgm_index_int32_find(pgm, 4, &value))
        printf("Key 4 found\n");
    else
        printf("Key 4 not found\n");

    dynamic_pgm_index_int32_find(pgm, 8, &value);
    printf("The value of key 8 is %d\n", value);

    printf("Range search [1, 10000) = ");
    dynamic_pgm_index_int32_iterator_t it = dynamic_pgm_index_int32_lower_bound(pgm, 1);
    while (dynamic_pgm_index_int32_iterator_next(pgm, it, &key, &value) && key < 10000)
        printf("(%d,%d), ", key, value);
    dynamic_pgm_index_int32_iterator_destroy(it);

    dynamic_pgm_index_int32_destroy(pgm);

    return 0;
}