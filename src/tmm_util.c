//
// Created by johnh on 8/8/2025.
//

#include <stdint.h>
#include <complex.h>
#include "tmm_util.h"


uint8_t tmm_array_zeros()
{

    return 0;
};


/**
 * @brief
 *
 * @param n Number of layers
 * @param mat nxn matrix from the outer scope
 * @return
 */
uint8_t tmm_matrix_zeros(const uint8_t n, double complex mat[n][n])
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mat[i][j] = 0.0 + 0.0 * I;
        }
    }
    return 0;
};


/**
 * @brief
 *
 * @param n Number of layers
 * @param mat nxn matrix from the outer scope
 * @param mat_copy nxn matrix from the outer scope
 * @return
 */
uint8_t tmm_matrix_copy(
    const uint8_t n, const double complex mat[n][n], double complex mat_copy[n][n]
)
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            mat_copy[i][j] = mat[i][j];
        }
    }

    return 0;
};
