#include <stdio.h>
#include <stdlib.h>
#include "vector.h"

double resize(vector* v, int capacity) {
    double* array = realloc(v -> array, sizeof(double) * capacity);
    if (array) v -> array = array, v -> capacity = capacity;
}

double append(vector* v, double item) {
    if (v -> capacity == v -> size)
        resize(v, v -> capacity * 2);
    v -> array[v -> size++] = item;
}

// int main() {
//     vector v;
//     v.capacity = 5, v.size = 0;
//     v.array = malloc(sizeof(double) * v.capacity);

//     for (int i = 0; i < 6; i++)
//         append(&v, i);
    
//     printf("Capacity: %i\nSize: %i\n", v.capacity, v.size);

//     for (int i = 0; i < 6; i++)
//         printf("%lf ", v.array[i]);
    
//     free(v.array);
// }
