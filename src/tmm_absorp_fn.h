//
// Created by johnh on 8/2/2025.
//

#ifndef TMM_ABSORP_FN_H
#define TMM_ABSORP_FN_H


#include "tmm_coherent.h"


// // Define a struct for a Point in 2D space
// struct Point {
//     int x;
//     int y;
// };
//
// int main() {
//     // Declare and initialize a Point
//     struct Point p1 = { .x = 10, .y = 20 };
//
//     // Access members
//     printf("Point coordinates: (%d, %d)\n", p1.x, p1.y);
//
//     return 0;
// }


/**
 * @brief Absorption in a given layer is a pretty simple analytical function.
 *
 * The sum of four exponentials.
 *
 *   a(z) = A1*exp(a1*z) + A2*exp(-a1*z)
 *          + A3*exp(1j*a3*z) + conj(A3)*exp(-1j*a3*z)
 *
 * where a(z) is absorption at depth z, with z=0 being the start of the layer,
 * and A1,A2,a1,a3 are real numbers, with a1>0, a3>0, and A3 is complex.
 * The class stores these five parameters, as well as d, the layer thickness.
 *
 * This gives absorption as a fraction of intensity coming towards the first
 * layer of the stack.
 *
 *
 * This struct holds integer values for the x and y coordinates.
 */
typedef struct {
    double d;
    double a1;
    double a2;
    double a3;
    double complex A1;
    double complex A2;
    double complex A3;
} AbsorpAnalyticFn;


uint8_t AbsorpAnalyticFn_create(AbsorpAnalyticFn self);


uint8_t AbsorpAnalyticFn_destroy(AbsorpAnalyticFn self);


uint8_t fill_in(
    AbsorpAnalyticFn* self, CohTmmData* coh_tmm_data, int layer
    // AbsorpAnalyticFn* self
);


uint8_t copy(const AbsorpAnalyticFn* self, AbsorpAnalyticFn* a);


uint8_t run(AbsorpAnalyticFn* self, double z, double complex* absorp);


uint8_t flip(AbsorpAnalyticFn* self);


uint8_t scale(AbsorpAnalyticFn* self, double factor);


uint8_t add(AbsorpAnalyticFn* self, const AbsorpAnalyticFn* b);


#endif //TMM_ABSORP_FN_H
