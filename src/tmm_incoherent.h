//
// Created by johnh on 8/9/2025.
//

#ifndef TMM_INCOHERENT_H
#define TMM_INCOHERENT_H


#include <stdint.h>
#include <complex.h>
#include "tmm_coherent.h"


/**
 * @brief Incoherent TMM Data struct
 *
 *
 */
typedef struct {
    double T;  // real transmissivity
    double R;  // real reflectivity
    double* VW_list;
    double* power_entering_list;
    double complex* stackFB_list;
    CohTmmData* coh_tmm_data_list;
    CohTmmData* coh_tmm_bdata_list;
    uint8_t num_inc_layers;
} IncTmmData;


#endif //TMM_INCOHERENT_H
