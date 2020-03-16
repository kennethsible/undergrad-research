#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lambertw.h"
#include "vector.h"

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
    // double x0 = 500; // Gaussian Pulse Center (Uniform)
    double x0 = 800; // Gaussian Pulse Center (Adaptive)
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
        psi[i] = (*f0)(u[i]), psi[i + Nu] = 0;

    double x, a; int n, k;
    for (n = 0; n < (Nv - 1); n++) {
        for (k = 0; k < (Nu - 1); k++) {
            x = ((v[n] - u[k + 1]) + (v[n + 1] - u[k]))/4;
            a = 1 - (hu * hv)/8 * V(r(x));
            psi[k + 1 + Nu] = a * (psi[k + Nu] + psi[k + 1]) - psi[k];
        }
        for (k = 0; k < Nu; k++)
            psi[k] = psi[k + Nu], psi[k + Nu] = 0;
    }

    FILE *file = fopen("psi.txt", "wb");
    // fwrite(psi, sizeof(double), sizeof(psi), file);
    for (k = 0; k < Nu; k++)
        fprintf(file, "%lf  %0.15f\n", u[k], psi[k]);
    fclose(file);

    free(u); free(v); free(psi);
}

static void adaptive(tuple* uspan, tuple* vspan, double hu, double (*f0)(double), double tol) {
    double ui = uspan -> a, uf = uspan -> b;
    double vi = vspan -> a, vf = vspan -> b;
    double* u = range(ui, uf, hu);
    vector v; v.capacity = 10, v.size = 0;
    v.array = malloc(sizeof(double) * v.capacity);
    append(&v, vi);
    int Nu = (uf - ui)/hu + 1, Nv = 1;
    double max_step = 100;

    double* psi   = malloc(sizeof(double) * Nu);
    double* psiHR = malloc(sizeof(double) * 2*Nu);
    double* psiLR = malloc(sizeof(double) * Nu);
    for (int i = 0; i < Nu; i++) {
        psi[i]   = psiLR[i] = 0;
        psiHR[i] = psiHR[i + Nu] = 0;
    } psi[0] = (*f0)(v.array[0]);

    int k, n = 0; double x, a, hv = 100;
    while (v.array[n] < vf) {
        psiLR[0] = psiHR[0] = psiHR[Nu] = (*f0)(v.array[n]);
        for (k = 0; k < (Nu - 1); k++) {
            x = ((v.array[n] - u[k + 1]) + ((v.array[n] + hv) - u[k]))/4;
            a = 1 - (hu * hv)/8 * V(r(x));
            psiLR[k + 1] = a * (psiLR[k] + psi[k + 1]) - psi[k];

            x = ((v.array[n] - u[k + 1]) + ((v.array[n] + hv/2) - u[k]))/4;
            a = 1 - (hu * hv/2)/8 * V(r(x));
            psiHR[k + 1] = a * (psiHR[k] + psi[k + 1]) - psi[k];

            x = (((v.array[n] + hv/2) - u[k + 1]) + ((v.array[n] + hv) - u[k]))/4;
            a = 1 - (hu * hv/2)/8 * V(r(x));
            psiHR[k + 1 + Nu] = a * (psiHR[k + Nu] + psiHR[k + 1]) - psiHR[k];
        }
        double epsilon = fabs(psiHR[Nu] - psiLR[0]);
        for (k = 1; k < Nu; k++) {
            double tmp = fabs(psiHR[k + Nu] - psiLR[k]);
            if (tmp > epsilon) epsilon = tmp;
        }
        if (epsilon < tol) {
            free(psi);
            psi = copy(psiLR, Nu);
            append(&v, v.array[n] + hv);
            n += 1; Nv += 1;
        }
        hv *= sqrt(tol/epsilon);
        if (hv > max_step) hv = max_step;
    }

    FILE *file = fopen("psi.txt", "wb");
    for (k = 0; k < Nu; k++)
        fprintf(file, "%lf  %0.15f\n", u[k], psi[k]);
    fclose(file);

    free(u); free(v.array); free(psi); free(psiHR); free(psiLR);
}

int main() {
    tuple uspan, vspan;
    uspan.a = 0, uspan.b = 1200;
    vspan.a = 0, vspan.b = 1600;
    double hu = 0.4;
    adaptive(&uspan, &vspan, hu, f0, 1e-7);
    return 0;
}
