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

