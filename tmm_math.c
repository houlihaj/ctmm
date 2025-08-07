//
// Created by johnh on 8/6/2025.
//

#include <stdio.h>
#include "tmm_math.h"


/**
 * @brief Makes a 2x2 array of [[a,b],[c,d]]
 *
 * Must declare the matrix in the outer scope and pass it into
 * this function for initialization of the elements
 *
 * @param mat1 2x2 array
 * @param mat2 2x2 array
 * @param product 2x2 array from the outer scope to hold the result
 * @return
 */
uint8_t tmm_matrix_product(
    const double complex mat1[2][2],
    const double complex mat2[2][2],
    double complex product[2][2]
)
{
    // Matrix multiplication product = mat1 * mat2
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            product[i][j] = 0.0 + 0.0 * I;
            for (int k = 0; k < 2; ++k) {
                product[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return 0;
}


/**
 * @brief Makes a 2x2 array of [[a,b],[c,d]]
 *
 * Must declare the matrix in the outer scope and pass it into
 * this function for initialization of the elements
 *
 * @param mat1 2x2 array from the outer scope; also the output result array
 * @param divisor the scalar
 * @return
 */
uint8_t tmm_scalar_division(
    double complex mat1[2][2], const double complex divisor
)
{
    if (divisor == 0.0) {
        fprintf(stderr, "Error: Division by zero is undefined.\n");
        return 1;  // return error
    }

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            mat1[i][j] /= divisor;
        }
    }
    return 0;
}


// #include <stdio.h>
// #include <complex.h>
//
// int main() {
//     // Define two 2x2 complex matrices A and B
//     double complex A[2][2] = {
//         {1.0 + 2.0*I, 3.0 + 4.0*I},
//         {5.0 + 6.0*I, 7.0 + 8.0*I}
//     };
//
//     double complex B[2][2] = {
//         {9.0 + 1.0*I, 2.0 + 3.0*I},
//         {4.0 + 5.0*I, 6.0 + 7.0*I}
//     };
//
//     // Result matrix C
//     double complex C[2][2];
//
//     // Matrix multiplication C = A * B
//     for (int i = 0; i < 2; ++i) {
//         for (int j = 0; j < 2; ++j) {
//             C[i][j] = 0.0 + 0.0*I;
//             for (int k = 0; k < 2; ++k) {
//                 C[i][j] += A[i][k] * B[k][j];
//             }
//         }
//     }
//
//     // Print the result
//     printf("Matrix C (result of A * B):\n");
//     for (int i = 0; i < 2; ++i) {
//         for (int j = 0; j < 2; ++j) {
//             printf("(%g%+gi) ", creal(C[i][j]), cimag(C[i][j]));
//         }
//         printf("\n");
//     }
//
//     return 0;
// }


