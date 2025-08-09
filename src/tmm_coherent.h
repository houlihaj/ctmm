//
// Created by jhoulihan on 8/5/2025.
//

#ifndef TMM_COHERENT_H
#define TMM_COHERENT_H

#include <stdint.h>
#include <complex.h>


/**
 * @brief asdfas
 *
 * This struct holds integer values for the x and y coordinates.
 */
typedef struct {
    double complex r;  // complex reflection amplitude (i.e. reflection coefficient)
    double complex t;  // complex transmission amplitude (i.e. transmission coefficient)
    double R;  // real reflectivity
    double T;  // real transmissivity
    uint8_t num_layers;
    double power_entering;
    double* vw_list;
    double complex* kz_list;
    double complex* th_list;
    double pol;
    double complex* n_list;
    double* d_list;
    double th_0;
    double lam_vac;
} CohTmmData;


uint8_t CohTmmData_create(CohTmmData* self, uint8_t num_layers);


uint8_t CohTmmData_destroy(CohTmmData* self);


#endif //TMM_COHERENT_H
