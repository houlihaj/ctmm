//
// Created by johnh on 8/2/2025.
//

#include <stdlib.h>
#include <stdint.h>
#include <complex.h>
#include <math.h>
#include "tmm_absorp_fn.h"
#include "tmm_coherent.h"


// /**
//  * @brief Absorption in a given layer is a pretty simple analytical function.
//  *
//  * The sum of four exponentials.
//  *
//  *   a(z) = A1*exp(a1*z) + A2*exp(-a1*z)
//  *          + A3*exp(1j*a3*z) + conj(A3)*exp(-1j*a3*z)
//  *
//  * where a(z) is absorption at depth z, with z=0 being the start of the layer,
//  * and A1,A2,a1,a3 are real numbers, with a1>0, a3>0, and A3 is complex.
//  * The class stores these five parameters, as well as d, the layer thickness.
//  *
//  * This gives absorption as a fraction of intensity coming towards the first
//  * layer of the stack.
//  *
//  *
//  * This struct holds integer values for the x and y coordinates.
//  */
// typedef struct {
//     double complex d;
//     double complex a1;
//     double complex a2;
//     double complex a3;
//     double complex A1;
//     double complex A2;
//     double complex A3;
// } AbsorpAnalyticFn;


/**
 * @brief
 *
 * @param self
 * @return
 */
uint8_t AbsorpAnalyticFn_create(AbsorpAnalyticFn self)
{
    return 0;
}


/**
 * @brief
 *
 * @param self
 * @return
 */
uint8_t AbsorpAnalyticFn_destroy(AbsorpAnalyticFn self)
{
    // free( (void*)self );  // TODO: must free the memory properly
    return 0;
}


/**
 * @brief
 *
 * @param
 * @param
 * @return
 */
uint8_t fill_in(
    // AbsorpAnalyticFn self, CohTmmData coh_tmm_data, int layer
    AbsorpAnalyticFn* self
)
{
    // Temporary!!!
    self->a1 = 1.0;
    self->d = 10.0;

    self->A1 = 1.0 + 0.0 * I;
    self->A2 = 1.2857 + 13.660 * I;  // Aluminum at 1.5 micron

    return 0;
}


/**
 * @brief Create copy of an AbsorpAnalyticFn object (i.e. C struct)
 *
 * @param self
 * @param a The AbsorpAnalyticFn copy
 * @return
 */
uint8_t copy(const AbsorpAnalyticFn self, AbsorpAnalyticFn a)
{
    a.A1 = self.A1;
    a.A2 = self.A2;
    a.A2 = self.A3;
    a.a1 = self.a1;
    a.a3 = self.a3;
    a.d = self.d;
    return 0;
}


/**
 * @brief Compute absorption at a given depth z, where z = 0 is the
 * start of the layer
 *
 * @param self
 * @param z
 * @return
 */
uint8_t run(AbsorpAnalyticFn self, const double z)
{
    return 0;
}


/**
 * @brief
 *
 * @param self
 * @return
 */
uint8_t flip(AbsorpAnalyticFn* self)
{
    const double complex newA1 = self->A2 * exp(-self->a1 * self->d);
    const double complex newA2 = self->A1 * exp(self->a1 * self->d);
    self->A1 = newA1;
    self->A2 = newA2;
    // self.A3 = conj(self.A3 * exp(1j * self.a3 * self.d));
    return 0;
}


/**
 * @brief Multiplies the absorption at each point by "factor".
 *
 * @param self
 * @param factor
 * @return
 */
uint8_t scale(AbsorpAnalyticFn self, const double factor)
{
    self.A1 *= factor;
    self.A2 *= factor;
    self.A3 *= factor;
    return 0;
}


/**
 * @brief Adds another compatible absorption analytical function
 *
 * @param self
 * @param b Another compatible absorption analytical function
 * @return
 */
uint8_t add(AbsorpAnalyticFn self, const AbsorpAnalyticFn b)
{
    self.A1 += b.A1;
    self.A2 += b.A2;
    self.A3 += b.A3;
    return 0;
}
