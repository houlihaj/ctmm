//
// Created by johnh on 8/9/2025.
//

#ifndef TMM_INCOHERENT_H
#define TMM_INCOHERENT_H


#include <stdint.h>
#include <complex.h>


/**
 * @brief asdfas
 *
 *
 */
typedef struct {
    double complex r;  // complex reflection amplitude (i.e. reflection coefficient)
    double complex t;  // complex transmission amplitude (i.e. transmission coefficient)
    double R;  // real reflectivity
    double T;  // real transmissivity
    double power_entering;
    double* vw_list;
    double complex* kz_list;
    double complex* th_list;
    double complex* n_list;
    double* d_list;
    double pol;
    double th_0;
    double lam_vac;
    uint8_t num_layers;
} IncTmmData;


#endif //TMM_INCOHERENT_H
