//
// Created by johnh on 8/9/2025.
//

#include <stdio.h>
#include <stdlib.h>
#include "tmm_incoherent.h"
#include "tmm_coherent.h"


/**
 * @brief Initialize IncTmmData structure passed by pointer.
 *
 * @param self
 * @param num_inc_layers
 * @return
 */
uint8_t IncTmmData_create(IncTmmData* self, const uint8_t num_inc_layers)
{
    self->num_inc_layers = num_inc_layers;

    self->VW_list = malloc((num_inc_layers * 2) * sizeof(double));
    self->power_entering_list = malloc(num_inc_layers * sizeof(double));
    self->stackFB_list = malloc(num_inc_layers * sizeof(double complex));
    self->coh_tmm_data_list = malloc(num_inc_layers * sizeof(CohTmmData));
    self->coh_tmm_bdata_list = malloc(num_inc_layers * sizeof(CohTmmData));

    if (
        !self->VW_list
        || !self->power_entering_list
        || !self->stackFB_list
        || !self->coh_tmm_data_list
        || !self->coh_tmm_bdata_list
    )
    {
        fprintf(stderr, "Error allocating memory for at least one array\n");
        free(self->VW_list);
        free(self->power_entering_list);
        free(self->stackFB_list);
        free(self->coh_tmm_data_list);
        free(self->coh_tmm_bdata_list);
        return 1;  // Return error code
    }

    return 0;
}


/**
 * @brief Free memory allocated for IncTmmData and its members.
 *
 * @param self
 * @return
 */
uint8_t IncTmmData_destroy(IncTmmData* self)
{
    free(self->VW_list);
    free(self->power_entering_list);
    free(self->stackFB_list);
    free(self->coh_tmm_data_list);
    free(self->coh_tmm_bdata_list);
    return 0;
}


/**
 * @brief Initialize IncTmmData structure passed by pointer.
 *
 * @param self
 * @param num_layers
 * @param num_inc_layers
 * @param num_stacks
 * @return
 */
uint8_t IncGroupLayersData_create(
    IncGroupLayersData* self,
    const uint8_t num_layers,
    const uint8_t num_inc_layers,
    const uint8_t num_stacks
)
{
    self->num_layers = num_layers;

    self->stack_d_list = malloc((num_stacks * 2) * sizeof(double));
    self->stack_n_list = malloc((num_stacks * 2) * sizeof(double complex));
    self->all_from_inc = malloc(num_inc_layers * sizeof(uint8_t));
    self->inc_from_all = malloc(num_layers * sizeof(uint8_t));
    self->all_from_stack = malloc(num_layers * sizeof(uint8_t));
    self->stack_from_all = malloc((num_layers * 2) * sizeof(uint8_t));
    self->inc_from_stack = malloc(num_stacks * sizeof(uint8_t));
    self->stack_from_inc = malloc(num_inc_layers * sizeof(uint8_t));

    if (
        !self->stack_d_list
        || !self->stack_n_list
        || !self->all_from_inc
        || !self->inc_from_all
        || !self->all_from_stack
        || !self->stack_from_all
        || !self->inc_from_stack
        || !self->stack_from_inc
    )
    {
        fprintf(stderr, "Error allocating memory for at least one array\n");
        free(self->stack_d_list);
        free(self->stack_n_list);
        free(self->all_from_inc);
        free(self->inc_from_all);
        free(self->all_from_stack);
        free(self->stack_from_all);
        free(self->inc_from_stack);
        free(self->stack_from_inc);
        return 1;  // Return error code
    }

    return 0;
}


/**
 * @brief Free memory allocated for IncTmmData and its members.
 *
 * @param self
 * @return
 */
uint8_t IncGroupLayersData_destroy(IncGroupLayersData* self)
{
    free(self->stack_d_list);
    free(self->stack_n_list);
    free(self->all_from_inc);
    free(self->inc_from_all);
    free(self->all_from_stack);
    free(self->stack_from_all);
    free(self->inc_from_stack);
    free(self->stack_from_inc);
    return 0;
}
