//
// Created by jhoulihan on 8/5/2025.
//

#include <stdio.h>
#include <stdlib.h>
#include "tmm_coherent.h"


/**
 * @brief Initialize CohTmmData structure passed by pointer.
 *
 * @param self
 * @param num_layers
 * @return
 */
uint8_t CohTmmData_create(CohTmmData* self, const uint8_t num_layers)
{
    self->num_layers = num_layers;

    self->vw_list = malloc(num_layers * sizeof(double));
    self->kz_list = malloc(num_layers * sizeof(double complex));
    self->th_list = malloc(num_layers * sizeof(double));
    self->n_list = malloc(num_layers * sizeof(double complex));
    self->d_list = malloc(num_layers * sizeof(double));

    if (!self->vw_list || !self->kz_list || !self->th_list || !self->n_list || !self->d_list) {
        fprintf(stderr, "Error allocating memory for at least one array\n");
        free(self->vw_list);
        free(self->kz_list);
        free(self->th_list);
        free(self->n_list);
        free(self->d_list);
        return 1;  // Return error code
    }

    return 0;
}


/**
 * @brief Free memory allocated for CohTmmData and its members.
 *
 * @param self
 * @return
 */
uint8_t CohTmmData_destroy(CohTmmData* self)
{
    free(self->vw_list);
    free(self->kz_list);
    free(self->th_list);
    free(self->n_list);
    free(self->d_list);
    return 0;
}
