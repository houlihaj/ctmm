//
// Created by johnh on 8/2/2025.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <complex.h>
#include <math.h>
#include "tmm_absorp_fn.h"
#include "tmm_coherent.h"


/**
 * @brief
 *
 * Done!!!
 *
 * @param
 * @param
 * @return
 */
uint8_t fill_in(
    AbsorpAnalyticFn* self, CohTmmData* coh_tmm_data, const int layer
)
{
    const uint8_t pol = coh_tmm_data->pol;
    // double v = coh_tmm_data->vw_list[layer][0];
    const double complex v = coh_tmm_data->vw_list[layer * 2];  // TODO: needs thorough validation
    // double w = coh_tmm_data->vw_list[layer][1];
    const double complex w = coh_tmm_data->vw_list[layer * 2 + 1];  // TODO: needs thorough validation
    const double complex kz = coh_tmm_data->kz_list[layer];
    const double complex n = coh_tmm_data->n_list[layer];
    const double complex n_0 = coh_tmm_data->n_list[0];
    const double th_0 = coh_tmm_data->th_0;
    const double complex th = coh_tmm_data->th_list[layer];
    self->d = coh_tmm_data->d_list[layer];

    self->a1 = 2 * cimag(kz);
    self->a3 = 2 * creal(kz);

    double temp;
    if (pol == 0)  // s-polarization
    {
        temp = (
            cimag(n * cos(th) * kz) / creal(n_0 * cos(th_0))
        );
        self->A1 = temp * creal(w) * creal(w) + cimag(w) * cimag(w);
        self->A2 = temp * creal(v) * creal(v) + cimag(v) * cimag(v);
        self->A3 = temp * v * conj(w);
    } else  // p-polarization
    {
        temp = (
            creal(2 * cimag(kz) * creal( n * cos(conj(th)) ))
            / creal( n_0 * conj(cos(th_0)) )
        );
        self->A1 = temp * creal(w) * creal(w) + cimag(w) * cimag(w);
        self->A2 = temp * creal(v) * creal(v) + cimag(v) * cimag(v);
        self->A3 = (
            v
            * conj(w)
            * (
                -2 * creal(kz) * cimag(n * cos(conj(th)) )
                / creal( n_0 * conj(cos(th_0)) )
            )
        );
    }

    return 0;
}


/**
 * @brief Create copy of an AbsorpAnalyticFn object (i.e. C struct)
 *
 * Done!!!
 *
 * @param self
 * @param a The AbsorpAnalyticFn copy
 * @return
 */
uint8_t copy(const AbsorpAnalyticFn* self, AbsorpAnalyticFn* a)
{
    a->A1 = self->A1;
    a->A2 = self->A2;
    a->A2 = self->A3;
    a->a1 = self->a1;
    a->a3 = self->a3;
    a->d = self->d;
    return 0;
}


/**
 * @brief Compute absorption at a given depth z, where z = 0 is the
 * start of the layer
 *
 * Done!!!
 *
 * @param self
 * @param z
 * @return
 */
uint8_t run(
    AbsorpAnalyticFn* self, const double z, double complex* absorp
)
{
    *absorp = (
        self->A1 * exp(self->a1 * z)
        + self->A2 * exp(-self->a1 * z)
        + self->A3 * cexp(1 * I * self->a3 * z)
        + conj(self->A3) * cexp(-1 * I * self->a3 * z)
    );
    return 0;
}


/**
 * @brief Flip the function front-to-back
 *
 * Done!!!
 *
 * To describe a(d - z) instead of a(z), where d is layer thickness.
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
    self->A3 = conj(self->A3 * cexp(I * self->a3 * self->d));
    return 0;
}


/**
 * @brief Multiplies the absorption at each point by "factor".
 *
 * Done!!!
 *
 * @param self
 * @param factor
 * @return
 */
uint8_t scale(AbsorpAnalyticFn* self, const double factor)
{
    self->A1 *= factor;
    self->A2 *= factor;
    self->A3 *= factor;
    return 0;
}


/**
 * @brief Adds another compatible absorption analytical function
 *
 * Done!!!
 *
 * @param self
 * @param b Another compatible absorption analytical function
 * @return
 */
uint8_t add(AbsorpAnalyticFn* self, const AbsorpAnalyticFn* b)
{
    if (b->a1 != self->a1 || b->a3 != self->a3)
    {
        printf("[ValueError] Incompatible absorption analytical functions!");
        return 1;  // error
    }
    self->A1 += b->A1;
    self->A2 += b->A2;
    self->A3 += b->A3;
    return 0;
}
