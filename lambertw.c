#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static double f(double x, double z) {
    return x * exp(x) - z;
}

static double df(double x) {
    return exp(x) * (x + 1);
}

/* lambertw(z, tol)
 *
 * Solves z = w * exp(w) for w using the Newton-Raphson method
 * with the specified tolerance, assuming the principle branch.
 */
double lambertw(double z) {
    double w = z < 1 ? z : log(z);
    if (f(w, z) == 0) return w;
    int iter = 0;
    while (fabs(f(w, z)) > 1e-8) {
        w += -f(w, z)/df(w);
        if ((iter++) > 1000) break;
    }
    return w;
}
