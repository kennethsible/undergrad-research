#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lambertw.h"
#include "vector.h"

#define psi(n, k) psi[n*Nu + k]

typedef struct {
    double a;
    double b;
} tuple;

static double* copy(double* array, int size) {
    double* duplicate = malloc(sizeof(double) * size);
    for (int i = 0; i < size; i++)
        duplicate[i] = array[i];
    return duplicate;
}

static double* range(double start, double stop, double step) {
    int size = (stop - start)/step + 1;
    double* x = malloc(sizeof(double) * size);
    if (x) { // Allocation Check
        for (int i = 0; i < size; i++)
            x[i] = start + i*step;
    } return x;
}

static double V(double r) { // Schwarzschild Potential
    int l = 3, r_s = 2; /* Angular Momentum, Schwarzschild Radius */
    return (1 - r_s/r)*(l*(l + 1)/(r*r) - 3*r_s/pow(r, 3));
}

static double r(double x) { // Tortoise Coordinate
    return 2*(lambertw(exp((x - 2)/2)) + 1);
}

static double f0(double x) {
    double x0 = 500; // Gaussian Pulse Center
    return 0.5*(exp(-0.1*(x - x0)*(x - x0)));
}

static void uniform(tuple* uspan, tuple* vspan, double hu, double hv, double (*f0)(double)) {
    double ui = uspan -> a, uf = uspan -> b;
    double vi = vspan -> a, vf = vspan -> b;
    double* u = range(ui, uf, hu);
    double* v = range(vi, vf, hv);
    int Nu = (uf - ui)/hu + 1, Nv = (vf - vi)/hv + 1;
    
    double* psi = malloc(sizeof(double) * 2*Nu);
    for (int i = 0; i < Nu; i++)
        psi(0, i) = (*f0)(u[i]), psi(1, i) = 0;

    double x, a; int n, k;
    for (n = 0; n < (Nv - 1); n++) {
        psi(1, 0) = psi(0, 0) + (hv/hu)*(psi(0, 1) - psi(0, 0));
        psi(1, Nu - 1) = psi(0, Nu - 1) - (hv/hu)*(psi(0, Nu - 1) - psi(0, Nu - 2));
        for (k = 0; k < (Nu - 2); k++) {
            x = ((v[n] - u[k + 1]) + (v[n + 1] - u[k]))/4;
            a = 1 - (hu * hv)/8 * V(r(x));
            psi(1, k + 1) = a * (psi(1, k) + psi(0, k + 1)) - psi(0, k);
        }
        for (k = 0; k < Nu; k++)
            psi(0, k) = psi(1, k), psi(1, k) = 0;
    }

    FILE *file = fopen("psi.txt", "wb");
    for (k = 0; k < Nu; k++)
        fprintf(file, "%lf %0.15f\n", u[k], psi(0, k));
    fclose(file);

    free(u); free(v); free(psi);
}

static void uniformHR(tuple* uspan, tuple* vspan, double hu, double hv, double (*f0)(double)) {
    double ui = uspan -> a, uf = uspan -> b;
    double vi = vspan -> a, vf = vspan -> b;
    double* u = range(ui, uf, hu);
    double* v = range(vi, vf, hv);
    int Nu = (uf - ui)/hu + 1, Nv = (vf - vi)/hv + 1;

    double* psi = malloc(sizeof(double) * 3*Nu);
    for (int i = 0; i < Nu; i++)
        psi(0, i) = psi(1, i) = (*f0)(u[i]), psi(2, i) = 0;
    double* v0 = range(vi, vi + hv, hv/100);

    double x, a; int n, k;
    for (n = 0; n < 100; n++) {
        psi(2, 0) = psi(1, 0) + (hv/hu)*(psi(1, 1) - psi(1, 0));
        psi(2, Nu - 1) = psi(1, Nu - 1) - (hv/hu)*(psi(1, Nu - 1) - psi(1, Nu - 2));
        for (k = 0; k < (Nu - 2); k++) {
            x = ((v0[n] - u[k + 1]) + (v0[n + 1] - u[k]))/4;
            a = 1 - (hu * hv)/8 * V(r(x));
            psi(2, k + 1) = a * (psi(2, k) + psi(1, k + 1)) - psi(1, k);
        }
        for (k = 0; k < Nu; k++)
            psi(1, k) = psi(2, k), psi(2, k) = 0;
    }

    for (n = 1; n < (Nv - 2); n++) {
        psi(2, 0) = psi(1, 0) + (hv/hu)*(psi(1, 1) - psi(1, 0));
        psi(2, 1) = psi(1, 1) + (hv/hu)*(psi(1, 2) - psi(1, 1));
        psi(2, Nu - 1) = psi(1, Nu - 1) - (hv/hu)*(psi(1, Nu - 1) - psi(1, Nu - 2));
        psi(2, Nu - 2) = psi(1, Nu - 2) - (hv/hu)*(psi(1, Nu - 2) - psi(1, Nu - 3));
        for (k = 0; k < (Nu - 4); k++) {
            x = ((v[n] - u[k + 1]) + (v[n + 1] - u[k]))/4;
            a = 1 - (hu * hv)/8 * V(r(x));
            psi(2, k + 2) = 8*(psi(1, k + 1) - psi(1, k) - psi(0, k + 1) + psi(0, k)) + (psi(2, k) + psi(0, k + 2) - psi(0, k)) + (hu * hv) * V(r(x)) * psi(0, k);
        }
        for (k = 0; k < Nu; k++) {
            psi(0, k) = psi(1, k);
            psi(1, k) = psi(2, k);
            psi(2, k) = 0;
        }
    }

    FILE *file = fopen("psi.txt", "wb");
    for (k = 0; k < Nu; k++)
        fprintf(file, "%lf %0.15f\n", u[k], psi(1, k));
    fclose(file);

    free(u); free(v); free(psi); free(v0);
}

int main() {
    tuple uspan, vspan;
    uspan.a = 400, uspan.b = 625;
    vspan.a = 0,   vspan.b = 7*0.4;
    double hu = 0.4, hv = 0.4;
    uniformHR(&uspan, &vspan, hu, hv, f0);

    // uspan.a = 400, uspan.b = 625;
    // vspan.a = 0,   vspan.b = 550;
    // double hu = 0.4, hv = 0.4;
    // uniform(&uspan, &vspan, hu, hv, f0);
    return 0;
}
