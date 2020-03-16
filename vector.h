#ifndef VECTOR_H
#define VECTOR_H

typedef struct vector {
    double* array;
    int capacity;
    int size;
} vector;

double resize(vector* v, int capacity);
double append(vector* v, double item);

#endif
