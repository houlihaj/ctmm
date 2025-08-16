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


/**
 * @brief Incoherent TMM Group Layers Data struct
 *
 *
 */
typedef struct {
    IncTmmData inc_tmm_data;  // how to properly allocate memory???
    double* stack_d_list;
    double complex* stack_n_list;
    uint8_t* all_from_inc;
    uint8_t* inc_from_all;
    uint8_t* all_from_stack;
    uint8_t* stack_from_all;
    uint8_t* inc_from_stack;
    uint8_t* stack_from_inc;
    uint8_t num_stacks;
    uint8_t num_inc_layers;
    uint8_t num_layers;
} IncGroupLayersData;


uint8_t IncTmmData_create(
    IncTmmData* self, uint8_t num_inc_layers
);


uint8_t IncTmmData_destroy(IncTmmData* self);


uint8_t IncGroupLayersData_create(
    IncGroupLayersData* self, uint8_t num_layers, uint8_t num_inc_layers, uint8_t num_stacks
);


uint8_t IncGroupLayersData_destroy(IncGroupLayersData* self);


#endif //TMM_INCOHERENT_H
