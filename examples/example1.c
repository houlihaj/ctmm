//
// Created by johnh on 8/2/2025.
//

#include <stdio.h>
#include <stdint.h>
#include <complex.h>
#include <math.h>

#include "src/tmm_absorp_fn.h"
#include "src/tmm_coherent.h"
#include "src/tmm_core.h"


int main(int argc, char** argv) {
    printf("Hello, World!\nMy name is John\n");

    // Compute Snell's law
    const double n_1 = 1.0;
    const double n_2 = 1.5;
    // const double th_1 = 0.17;  // units of rad
    const double th_1 = 0.0;  // units of rad
    double th_2_guess;  // units of rad
    snell(n_1, n_2, th_1, &th_2_guess);
    printf("th_2_guess (rad): %lf\n", th_2_guess);

    // Compute the complex reflection and transmission coefficients
    // (i.e. reflection and transmission amplitudes)
    const double complex eta_1 = 1.0 + 0.0 * I;
    const double complex eta_2 = 1.2857 + 13.660 * I;  // Aluminum at 1.5 micron
    double complex r;
    interface_r(0, eta_1, eta_2, th_1, th_2_guess, &r);
    printf("reflection coefficient, r = %.3f + %.3fi\n", creal(r), cimag(r));
    double complex t;
    interface_t(0, eta_1, eta_2, th_1, th_2_guess, &t);
    printf("transmission coefficient, t = %.3f + %.3fi\n", creal(t), cimag(t));

    // Compute reflectivity R and transmittivity
    double R;
    R_from_r(r, &R);
    printf("reflectivity, R = %.3f\n", R);

    double T;
    T_from_t(0, t, eta_1, eta_2, th_1, th_2_guess, &T);
    printf("transmittivity, R = %.3f\n", T);

    ////////////////////////////////////////////////////////////////////////////
    // // Test AbsorpAnalyticFn
    //
    // AbsorpAnalyticFn absorp_fn;
    // fill_in(&absorp_fn);
    // printf(
    //     "absorp_fn.A1, absorp_fn.A2: (%.3f + %.3fi, %.3f + %.3fi)\n",
    //     creal(absorp_fn.A1), cimag(absorp_fn.A1), creal(absorp_fn.A2), cimag(absorp_fn.A2)
    // );
    // flip(&absorp_fn);
    // printf(
    //     "absorp_fn.A1, absorp_fn.A2: (%.3f + %.3fi, %.3f + %.3fi)\n",
    //     creal(absorp_fn.A1), cimag(absorp_fn.A1), creal(absorp_fn.A2), cimag(absorp_fn.A2)
    // );

    ////////////////////////////////////////////////////////////////////////////
    // // Test CohTmmData instantiation, memory allocation, and deallocation
    //
    // CohTmmData coh_tmm_data;  // allocate struct in outer scope
    // const uint8_t num_layers = 3;
    // CohTmmData_create(&coh_tmm_data, num_layers);
    //
    // double complex n_list[] = {  // number of items must equal num_layers
    //     1.0 + 0.0 * I, 2.5 + 0.0 * I, 1.5 + 0.0 * I
    // };
    // for (int i = 0; i < num_layers; i++)
    // {
    //     coh_tmm_data.n_list[i] = n_list[i];
    // }
    // for (int i = 0; i < num_layers; i++)
    // {
    //     printf(
    //         "complex index, n = %.3f + %.3fi\n",
    //         creal(coh_tmm_data.n_list[i]), cimag(coh_tmm_data.n_list[i])
    //     );
    // }
    //
    // coh_tmm_data.r = r;  // set reflection amplitude
    // printf("reflection coefficient, r = %.3f + %.3fi\n", creal(coh_tmm_data.r), cimag(coh_tmm_data.r));
    //
    // coh_tmm_data.R = R;  // set reflectivity
    // printf("reflectivity, R = %.3f\n", coh_tmm_data.R);
    //
    // CohTmmData_destroy(&coh_tmm_data);

    ////////////////////////////////////////////////////////////////////////////
    // Execute coh_tmm() to compute reflection and transmission coefficients

    CohTmmData coh_tmm_data;  // allocate struct in outer scope
    const uint8_t num_layers = 4;
    CohTmmData_create(&coh_tmm_data, num_layers);

    // list of layer thicknesses in nm
    const double d_list[] = {INFINITY, 100, 300, INFINITY};  // must start and end with INFINITY
    // list of refractive indices
    const double complex n_list[] = {
        1.0 + 0.0 * I, 2.2 + 0.0 * I, 3.3 + 0.3 * I, 1.0 + 0.0 * I
    };

    coh_tmm(0, n_list, d_list, num_layers, 0.0, 1550.0, &coh_tmm_data);

    printf("reflection coefficient, r = %.3f + %.3fi\n", creal(coh_tmm_data.r), cimag(coh_tmm_data.r));
    printf("reflectivity, R = %.3f\n", coh_tmm_data.R);

    printf("transmission coefficient, t = %.3f + %.3fi\n", creal(coh_tmm_data.t), cimag(coh_tmm_data.t));
    printf("transmissivity, T = %.3f\n", coh_tmm_data.T);

    CohTmmData_destroy(&coh_tmm_data);

    return 0;
}
