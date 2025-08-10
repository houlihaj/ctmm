//
// Created by johnh on 8/8/2025.
//

#ifndef TMM_UTIL_H
#define TMM_UTIL_H


uint8_t tmm_make_2x2_array(
    double complex a,
    double complex b,
    double complex c,
    double complex d,
    double complex matrix[2][2]
);


uint8_t tmm_array_zeros();


uint8_t tmm_matrix_zeros(uint8_t n, double complex mat[n][n]);


uint8_t tmm_matrix_copy(
    uint8_t n, const double complex mat[n][n], double complex mat_copy[n][n]
);


#endif //TMM_UTIL_H
