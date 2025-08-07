//
// Created by johnh on 8/6/2025.
//

#ifndef TMM_MATH_H
#define TMM_MATH_H

#include <stdint.h>
#include <complex.h>


uint8_t tmm_matrix_product(
    const double complex mat1[2][2],
    const double complex mat2[2][2],
    double complex product[2][2]
);


uint8_t tmm_scalar_division(
    double complex mat1[2][2], const double complex divisor
);


#endif //TMM_MATH_H
