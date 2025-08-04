//
// Created by johnh on 8/2/2025.
//

#include <stdio.h>
#include <stdint.h>
#include <complex.h>
#include "tmm_core.h"


int main(int argc, char** argv) {
    printf("Hello, World!\nMy name is John\n");

    // Compute Snell's law
    const double n_1 = 1.0;
    const double n_2 = 1.5;
    const double th_1 = 0.17;  // units of rad
    double th_2_guess;  // units of rad
    snell(n_1, n_2, th_1, &th_2_guess);
    printf("th_2_guess (rad): %lf\n", th_2_guess);

    // Compute the complex reflection coefficient (i.e. reflection amplitude)
    const double complex eta1 = 1.0 + 0.0 * I;
    const double complex eta2 = 1.2857 + 13.660 * I;  // Aluminum at 1.5 micron
    double complex r;
    interface_r(0, eta1, eta2, th_1, th_2_guess, &r);
    printf("reflection coefficient, r = %.3f + %.3fi\n", creal(r), cimag(r));

    // Compute reflectivity

    return 0;
}
