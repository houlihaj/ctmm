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
#include "src/tmm_incoherent.h"


int main(int argc, char** argv) {
    printf("Hello, World!\nMy name is John\n");

    // Incoherent group layers
    const uint8_t num_layers = 7;
    double complex n_list[] = {
        1.0 + 0.0 * I, 1.5 + 1.0 * I, 1.5 + 2.0 * I, 2.5 + 1.0 * I, 1.5 + 1.5 * I, 3.5 + 1.0 * I, 1.0 + 0.0 * I
    };
    double d_list[] = {
        INFINITY, 100, 300, 200, 500, 100, INFINITY
    };
    uint8_t c_list[] = {
        0, 1, 1, 0, 0, 1, 0
    };
    const uint8_t num_inc_layers = 4;
    const uint8_t num_stacks = 2;

    IncGroupLayersData inc_group_layers_data;
    IncGroupLayersData_create(
        &inc_group_layers_data, num_layers, num_inc_layers, num_stacks
    );

    inc_group_layers(
        &inc_group_layers_data, n_list, d_list, c_list, num_layers, num_inc_layers, num_stacks
    );

    IncGroupLayersData_destroy(&inc_group_layers_data);

    return 0;
}
