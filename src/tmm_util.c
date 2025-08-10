//
// Created by johnh on 8/8/2025.
//

#include <stdint.h>
#include <complex.h>
#include "tmm_util.h"


/**
 * @brief Makes a 2x2 array of [[a,b],[c,d]]
 *
 * Done!!!
 *
 * Must declare the matrix in the outer scope and pass it into
 * this function for initialization of the elements
 *
 * @param a
 * @param b
 * @param c
 * @param d
 * @param matrix 2x2 array
 * @return
 */
uint8_t tmm_make_2x2_array(
    const double complex a,
    const double complex b,
    const double complex c,
    const double complex d,
    double complex matrix[2][2]
)
{
    matrix[0][0] = a;
    matrix[0][1] = b;
    matrix[1][0] = c;
    matrix[1][1] = d;

    return 0;
}


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
