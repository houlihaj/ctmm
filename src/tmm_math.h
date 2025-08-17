//
// Created by johnh on 8/6/2025.
//

#ifndef TMM_MATH_H
#define TMM_MATH_H

#include <stdint.h>
#include <complex.h>


#ifdef __cplusplus
extern "C" {
#endif


uint8_t tmm_matrix_product(
    const double complex mat1[2][2],
    const double complex mat2[2][2],
    double complex product[2][2]
);


uint8_t tmm_matrix_by_vector(
    const double complex mat1[2][2],
    const double complex mat2[2][1],  // copilot suggestion
    double complex product[2][1]  // copilot suggestion
);


uint8_t tmm_scalar_division(
    double complex mat1[2][2], const double complex divisor
);


uint8_t tmm_scalar_product(
    double complex mat1[2][2], const double complex scalar
);


uint8_t tmm_transpose(
    double complex vw[2][1], double complex vw_tr[1][2]
);


#ifdef __cplusplus
}
#endif


#endif //TMM_MATH_H
