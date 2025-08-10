//
// Created by johnh on 8/2/2025.
//

/**
 * Definitions:
 *
 *
 */

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <complex.h>
#include <float.h>
#include <math.h>
#include "tmm_absorp_fn.h"
#include "tmm_coherent.h"
#include "tmm_core.h"
#include "tmm_math.h"
#include "tmm_util.h"


typedef enum
{
    SUCCESS,
    ERROR,
    NOT_IMPLEMENTED
} TmmStatus;


/**
 * @brief Compute whether-or-not this is the forward-traveling wave
 *
 * Done!!!
 *
 * @param n Medium with complex refractive index, n
 * @param theta Angle of incidence
 * @param answer Pass by reference from outer scope
 * @return
 */
uint8_t is_forward_angle(
    const double complex n, const double theta, bool* answer
)
{
    const double complex ncostheta = n * cos(theta);

    if ( fabs( cimag(ncostheta) ) > 100.0 * DBL_EPSILON )
    {
        *answer = cimag(ncostheta) > 0;
    } else
    {
        *answer = creal(ncostheta) > 0;
    }

    // double-check the answer ... cannot be too careful!

    if ( answer )
    {

        if ( !( cimag(ncostheta) > -100.0 * DBL_EPSILON) )
        {
            printf(
                "It is not clear which beam is incoming vs outgoing. "
                "Weird index maybe?\n"
                "n: %.2f + %.2fi   angle: %.4f\n", creal(n), cimag(n), theta
            );
            return 1;
        }

        if ( !( creal(ncostheta) > -100.0 * DBL_EPSILON) )
        {
            printf(
                "It is not clear which beam is incoming vs outgoing. "
                "Weird index maybe?\n"
                "n: %.2f + %.2fi   angle: %.4f\n", creal(n), cimag(n), theta
            );
            return 1;
        }

        if ( !( creal( n * cos(conj(theta)) ) > -100.0 * DBL_EPSILON ) )
        {
            printf(
                "It is not clear which beam is incoming vs outgoing. "
                "Weird index maybe?\n"
                "n: %.2f + %.2fi   angle: %.4f\n", creal(n), cimag(n), theta
            );
            return 1;
        }

    } else
    {

        if ( !( cimag(ncostheta) < 100.0 * DBL_EPSILON) )
        {
            printf(
                "It is not clear which beam is incoming vs outgoing. "
                "Weird index maybe?\n"
                "n: %.2f + %.2fi   angle: %.4f\n", creal(n), cimag(n), theta
            );
            return 1;
        }

        if ( !( creal(ncostheta) < 100.0 * DBL_EPSILON) )
        {
            printf(
                "It is not clear which beam is incoming vs outgoing. "
                "Weird index maybe?\n"
                "n: %.2f + %.2fi   angle: %.4f\n", creal(n), cimag(n), theta
            );
            return 1;
        }

        if ( !( creal( n * cos(conj(theta)) ) < 100.0 * DBL_EPSILON ) )
        {
            printf(
                "It is not clear which beam is incoming vs outgoing. "
                "Weird index maybe?\n"
                "n: %.2f + %.2fi   angle: %.4f\n", creal(n), cimag(n), theta
            );
            return 1;
        }

    }

    return 0;
}


/**
 * @brief Return the angle theta in layer 2 with refractive index n_2
 *
 * Done!!!
 *
 * @param n_1
 * @param n_2
 * @param th_1
 * @param th_2_guess Pass by reference from outer scope
 * @return
 */
uint8_t snell(
    const double n_1,
    const double complex n_2,
    const double th_1,
    double* th_2_guess
)
{
    bool answer;
    is_forward_angle(n_2, *th_2_guess, &answer);  // must de-reference th_2_guess

    if ( answer )
    {
        *th_2_guess = asin(n_1 * sin(th_1) / n_2);  // units of rad
    } else
    {
        *th_2_guess = PI - asin(n_1 * sin(th_1) / n_2);  // units of rad
    }

    return 0;
}


/**
 * @brief
 *
 * Done!!!
 *
 * @param n_list
 * @param n_list_size
 * @param th_0
 * @param angles
 * @return
 */
uint8_t list_snell(
    const double complex n_list[],
    const uint8_t n_list_size,
    const double th_0,
    double complex angles[]
)
{
    // Important that the arcsin here is numpy.lib.scimath.arcsin, not
    // numpy.arcsin! (They give different results e.g. for arcsin(2).)
    const double n_list_0 = n_list[0];
    for (uint8_t i = 0; i < n_list_size; i++)
    {
        // TODO: May have to do creal(n_list_0) instead of n_list_0
        angles[i] = asin( n_list_0 * sin(th_0) / n_list[i] );
        // angles[i] = asin( creal(n_list_0) * sin(th_0) / n_list[i] );
    }

    // The first and last entry need to be the forward angle (the intermediate
    // layers do not matter, see https://arxiv.org/abs/1603.02720 Section 5)
    bool answer;
    is_forward_angle(n_list[0], angles[0], &answer);
    if (!answer)
    {
        double complex angle = angles[0];
        angles[0] = PI - angle;
    }

    is_forward_angle(n_list[n_list_size - 1], angles[n_list_size - 1], &answer);
    if (!answer)
    {
        double complex angle = angles[n_list_size - 1];
        angles[n_list_size - 1] = PI - angle;
    }

    return 0;
}


/**
 * @brief Complex reflection amplitude from Fresnel equations
 *
 * Done!!!
 *
 * @param polarization
 * @param n_i
 * @param n_f
 * @param th_i
 * @param th_f
 * @param r Pass by reference from outer scope
 * @return
 */
uint8_t interface_r(
    const uint8_t polarization,
    const double complex n_i,
    const double complex n_f,
    const double complex th_i,
    const double complex th_f,
    double complex* r
)
{
    if (polarization == 0)  // s-polarization
    {
        *r = (
            ( n_i * cos(th_i) - n_f * cos(th_f) )
            / ( n_i * cos(th_i) + n_f * cos(th_f) )
        );
    } else  // p-polarization
    {
        *r = (
            ( n_f * cos(th_i) - n_i * cos(th_f) )
            / ( n_f * cos(th_i) + n_i * cos(th_f) )
        );
    }

    return 0;
}


/**
 * @brief
 *
 * Done!!!
 *
 * @param polarization
 * @param n_i
 * @param n_f
 * @param th_i
 * @param th_f
 * @param t
 * @return
 */
uint8_t interface_t(
    const uint8_t polarization,
    const double complex n_i,
    const double complex n_f,
    const double complex th_i,
    const double complex th_f,
    double complex* t
)
{
    if (polarization == 0)  // s-polarization
    {
        *t = (
            2 * n_i * cos(th_i)
            / ( n_i * cos(th_i) + n_f * cos(th_f) )
        );
    } else  // p-polarization
    {
        *t = (
            2 * n_i * cos(th_i)
            / ( n_f * cos(th_i) + n_i * cos(th_f) )
        );
    }

    return 0;
}


/**
 * @brief Calculate reflected power R, starting with reflection amplitude r.
 *
 * Done!!!
 *
 * @param r
 * @param R
 * @return
 */
uint8_t R_from_r(const double complex r, double* R)
{
    // Method 1: Manually compute the modulus squared
    const double real_part = creal(r);
    const double imag_part = cimag(r);
    *R = real_part * real_part + imag_part * imag_part;

    return 0;
}


/**
 * @brief Calculate transmitted power T, starting with transmission amplitude t.
 *
 * Done!!!
 *
 * @param pol
 * @param t
 * @param n_i
 * @param n_f
 * @param th_i
 * @param th_f
 * @param T Pass by reference from outer scope
 * @return
 */
uint8_t T_from_t(
    const uint8_t pol,
    const double complex t,
    const double complex n_i,
    const double complex n_f,
    const double complex th_i,
    const double complex th_f,
    double* T
)
{
    // Method 1: Manually compute the modulus squared
    const double real_part = creal(t);
    const double imag_part = cimag(t);

    if (pol == 0)  // s-polarization
    {
        *T = (
            (real_part * real_part + imag_part * imag_part)
            * ( creal( n_f * cos(th_f) ) / creal( n_i * cos(th_i) ) )
        );
    } else  // p-polarization
    {
        *T = (
            (real_part * real_part + imag_part * imag_part)
            * ( creal( n_f * conj(cos(th_f)) ) / creal( n_i * conj(cos(th_i)) ) )
        );
    }

    return 0;
}


/**
 * @brief Calculate the power entering the first interface, starting with
 * complex reflection amplitude, r.
 *
 * Done!!!
 *
 * @param pol
 * @param r
 * @param n_i The refractive index of incident medium
 * @param th_i The complex propagation angle through incident medium
 * @param power Pass by reference from outer scope
 * @return
 */
uint8_t power_entering_from_r(
    const uint8_t pol,
    const double complex r,
    const double complex n_i,
    const double th_i,
    double* power
)
{
    if (pol == 0)  // s-polarization
    {
        *power = (
            creal( n_i * cos(th_i) * (1 + conj(r)) * (1 - r) )
            / creal( n_i * cos(th_i) )
        );
    } else  // p-polarization
    {
        *power = (
            creal( n_i * conj(cos(th_i)) * (1 + r) * (1 - conj(r)) )
            / creal( n_i * conj(cos(th_i)) )
        );
    }

    return 0;
}


/**
 * @brief Fraction of light intensity reflected at an interface
 *
 * Done!!!
 *
 * @param
 * @param
 * @return
 */
uint8_t interface_R(
    const uint8_t polarization,
    const double complex n_i,
    const double complex n_f,
    const double th_i,
    const double th_f,
    double* R
)
{
    double complex r;
    interface_r(polarization, n_i, n_f, th_i, th_f, &r);
    R_from_r(r, R);
    return 0;
}


/**
 * @brief Fraction of light intensity transmitted at an interface
 *
 * Done!!!
 *
 * @param
 * @param
 * @return
 */
uint8_t interface_T(
    const uint8_t polarization,
    const double complex n_i,
    const double complex n_f,
    const double th_i,
    const double th_f,
    double* T
)
{
    double complex t;
    interface_t(polarization, n_i, n_f, th_i, th_f, &t);
    T_from_t(polarization, t, n_i, n_f, th_i, th_f, T);
    return 0;
}


/**
 * @brief Main "coherent transfer matrix method" computation.
 *
 * Done!!!
 *
 * @param pol Light polarization, "s" or "p"
 * @param n_list
 * @param d_list
 * @param num_layers
 * @param th_0 The angle of incidence
 * @param lam_vac The vacuum wavelength of the light
 * @param coh_tmm_data Pass by reference from outer scope
 * @return
 */
uint8_t coh_tmm(
    const uint8_t pol,
    const double complex n_list[],
    const double d_list[],
    const uint8_t num_layers,
    const double th_0,
    const double lam_vac,
    CohTmmData* coh_tmm_data
)
{
    // d_list must start and end with INFINITY
    if (d_list[0] != INFINITY || d_list[num_layers - 1] != INFINITY)
    {
        printf("d_list must start and end with inf!\n");
        // return with error
        return 1;
    }

    // Must be a forward angle
    bool answer;
    is_forward_angle(n_list[0], th_0, &answer);
    if (!answer)
    {
        printf("Error in n0 or th0!\n");
        return 1;
    }

    // th_list is a list with, for each layer, the angle that the light
    // travels through the layer. Computed with Snell's law. Note that
    // the "angles" may be complex!
    double complex th_list[num_layers];
    list_snell(n_list, num_layers, th_0, th_list);

    // kz is the z-component of (complex) angular wavevector for
    // forward-moving wave. Positive imaginary part means decaying.
    double complex kz_list[num_layers];
    for (int i = 0; i < num_layers; i++)
    {
        kz_list[i] = (
            2 * M_PI * n_list[i] * cos(th_list[i]) / lam_vac
        );
    }

    // delta is the total phase accrued by traveling through a given layer.
    // Ignore warning about inf multiplication
    double complex delta[num_layers];
    for (int i = 0; i < num_layers; i++)
    {
        // TODO: when d_list[i] is INFINITY, may have to take creal(kz_list[i])
        delta[i] = kz_list[i] * d_list[i];
    }

    // For a very opaque layer, reset delta to avoid divide-by-0 and similar
    // errors. The criterion imag(delta) > 35 corresponds to single-pass
    // transmission < 1e-30 --- small enough that the exact value does not
    // matter.
    for (int i = 1; i < num_layers - 1; i++)
    {
        if ( cimag(delta[i]) > 35 )
        {
            const double complex delta_i = delta[i];
            delta[i] = creal(delta_i) + cimag(35);

            // TODO: Must add opacity check and warning!!!
            // TODO: function can be validated before implementing check
        }
    }

    // t_list[i,j] and r_list[i,j] are transmission and reflection amplitudes,
    // respectively, coming from i, going to j. Only need to calculate this when
    // j=i+1. (2D array is overkill but helps avoid confusion.)
    double complex t_list[num_layers][num_layers];
    tmm_matrix_zeros(num_layers, t_list);
    double complex r_list[num_layers][num_layers];
    tmm_matrix_zeros(num_layers, r_list);
    for (int i = 0; i < num_layers - 1; i++)
    {
        double complex t;
        interface_t(
            pol, n_list[i], n_list[i + 1], th_list[i], th_list[i + 1], &t
        );
        t_list[i][i + 1] = t;

        double complex r;
        interface_r(
            pol, n_list[i], n_list[i + 1], th_list[i], th_list[i + 1], &r
        );
        r_list[i][i + 1] = r;
    }

    // At the interface between the (n-1)st and nth material, let v_n be the
    // amplitude of the wave on the nth side heading forwards (away from the
    // boundary), and let w_n be the amplitude on the nth side heading backwards
    // (towards the boundary). Then (v_n,w_n) = M_n (v_{n+1},w_{n+1}). M_n is
    // M_list[n]. M_0 and M_{num_layers-1} are not defined.
    // My M is a bit different than Sernelius's, but Mtilde is the same.
    double complex M_list[num_layers][2][2];
    for (int i = 0; i < num_layers; i++)
    {
        tmm_matrix_zeros(2, M_list[i]);
    }
    for (int i = 1; i < num_layers - 1; i++)
    {
        double complex mat1[2][2];
        tmm_make_2x2_array( cexp(-1 * I * delta[i]), 0, 0, cexp(1 * I * delta[i]), mat1 );

        double complex mat2[2][2];
        tmm_make_2x2_array( 1, r_list[i][i + 1], r_list[i][i + 1], 1, mat2 );

        tmm_matrix_product(mat1, mat2, M_list[i]);
        tmm_scalar_product(M_list[i], (1 / t_list[i][i + 1]));
        // tmm_scalar_division(M_list[i], t_list[i][i + 1]);  // same results as scalar product
    }
    double complex Mtilde[2][2];
    tmm_make_2x2_array(1, 0, 0, 1, Mtilde);
    double complex Mtilde_copy[2][2];
    for (int i = 1; i < num_layers - 1; i++)
    {
        tmm_matrix_copy(2, Mtilde, Mtilde_copy);
        tmm_matrix_product(Mtilde_copy, M_list[i], Mtilde);
    }
    double complex mat1[2][2];
    tmm_make_2x2_array( 1, r_list[0][1], r_list[0][1], 1, mat1 );
    tmm_scalar_division(mat1, t_list[0][1]);
    tmm_matrix_copy(2, Mtilde, Mtilde_copy);
    tmm_matrix_product(mat1, Mtilde_copy, Mtilde);

    // Net complex transmission and reflection amplitudes
    const double complex r_coeff = Mtilde[1][0] / Mtilde[0][0];
    const double complex t_coeff = 1 / Mtilde[0][0];

    // vw_list[n] = [v_n, w_n]. v_0 and w_0 are undefined because
    // the 0th medium has no left interface.
    double complex vw_list[num_layers][2];
    // 2x1 array (2 rows, 1 column)
    double complex vw[2][1] = {{t_coeff},{0.0 + 0.0 * I}};
    // 1x2 array (1 rows, 2 column)
    double complex vw_tr[1][2];
    tmm_transpose(vw, vw_tr);
    vw_list[num_layers - 1][0] = vw_tr[0][0];
    vw_list[num_layers - 1][1] = vw_tr[0][1];

    for (int i = num_layers - 2; i > 0; i--)
    {
        tmm_matrix_by_vector(M_list[i], vw, vw);
        tmm_transpose(vw, vw_tr);
        vw_list[i][0] = vw_tr[0][0];
        vw_list[i][1] = vw_tr[0][1];
    }

    // Net transmitted and reflected power, as a proportion of the
    // incoming light power.
    double R;
    R_from_r(r_coeff, &R);
    double T;
    T_from_t(
        pol, t_coeff, n_list[0], n_list[num_layers - 1], th_list[0], th_list[num_layers - 1], &T
    );
    double power_entering;
    power_entering_from_r(pol, r_coeff, n_list[0], th_0, &power_entering);

    // Store the data in the struct
    for (int i = 0; i < num_layers; i++)
    {
        // coh_tmm_data->vw_list[i] = vw_list[i];  // vw_list is an array whose items are an array of size 2
        coh_tmm_data->kz_list[i] = kz_list[i];
        coh_tmm_data->th_list[i] = th_list[i];
        coh_tmm_data->n_list[i] = n_list[i];  // assign items in n_list to coh_tmm_data->n_list
        coh_tmm_data->d_list[i] = d_list[i];
    }
    coh_tmm_data->r = r_coeff;
    coh_tmm_data->t = t_coeff;
    coh_tmm_data->R = R;
    coh_tmm_data->T = T;
    coh_tmm_data->power_entering = power_entering;
    coh_tmm_data->pol = pol;
    coh_tmm_data->th_0 = th_0;
    coh_tmm_data->lam_vac = lam_vac;

    return 0;
}


/**
 * @brief Reverses the order of the stack and then runs coh_tmm()
 *
 * Done!!!
 *
 * @param pol
 * @param n_list
 * @param d_list
 * @param num_layers
 * @param th_0
 * @param lam_vac
 * @param coh_tmm_data Pass by reference from outer scope
 * @return
 */
uint8_t coh_tmm_reverse(
    const uint8_t pol,
    const double complex n_list[],
    const double d_list[],
    const uint8_t num_layers,
    const double th_0,
    const double lam_vac,
    CohTmmData* coh_tmm_data
)
{
    double th_f;
    snell(n_list[0], n_list[num_layers - 1], th_0, &th_f);

    double complex n_list_rev[num_layers];
    double d_list_rev[num_layers];
    for (int i = 0; i < num_layers; i++)
    {
        n_list_rev[i] = n_list[(num_layers - 1) - i];
        d_list_rev[i] = d_list[(num_layers - 1) - i];
    }

    coh_tmm(pol, n_list_rev, d_list_rev, num_layers, th_f, lam_vac, coh_tmm_data);

    return 0;
}


/**
 * @brief Calculates ellipsometric parameters, in radians.
 *
 * Done!!!
 *
 * Warning: Conventions differ. You may need to subtract pi/2 or whatever.
 *
 * @param n_list
 * @param d_list
 * @param num_layers
 * @param th_0
 * @param lam_vac
 * @param ellips_data Pass by reference from outer scope
 * @return
 */
uint8_t ellips(
    const double complex n_list[],
    const double d_list[],
    const uint8_t num_layers,
    const double th_0,
    const double lam_vac,
    EllipsData* ellips_data
)
{
    CohTmmData coh_tmm_data_s;
    coh_tmm(0, n_list, d_list, num_layers, th_0, lam_vac, &coh_tmm_data_s);
    const double complex rs = coh_tmm_data_s.r;

    CohTmmData coh_tmm_data_p;
    coh_tmm(1, n_list, d_list, num_layers, th_0, lam_vac, &coh_tmm_data_p);
    const double complex rp = coh_tmm_data_p.r;

    ellips_data->psi = (
        atan( creal(rp / rs) * creal(rp / rs) + cimag(rp / rs) * cimag(rp / rs) )
    );
    ellips_data->delta = carg(-rp / rs);  // carg to compute argument of a complex number
    return 0;
}


/**
 * @brief Calculates reflected and transmitted power for unpolarized light.
 *
 * Done!!!
 *
 * @param n_list
 * @param d_list
 * @param num_layers
 * @param th_0
 * @param lam_vac
 * @param R
 * @param T
 * @return
 */
uint8_t unpolarized_RT(
    const double complex n_list[],
    const double d_list[],
    const uint8_t num_layers,
    const double th_0,
    const double lam_vac,
    double* R,
    double* T
)
{
    CohTmmData coh_tmm_data_s;
    coh_tmm(0, n_list, d_list, num_layers, th_0, lam_vac, &coh_tmm_data_s);
    const double Rs = coh_tmm_data_s.R;
    const double Ts = coh_tmm_data_s.T;

    CohTmmData coh_tmm_data_p;
    coh_tmm(1, n_list, d_list, num_layers, th_0, lam_vac, &coh_tmm_data_p);
    const double Rp = coh_tmm_data_p.R;
    const double Tp = coh_tmm_data_p.T;

    *R = (Rs + Rp) / 2;
    *T = (Ts + Tp) / 2;

    return 0;
}


/**
 * @brief Calculate the Poynting vector, absorbed energy density, and E-field at a specific location.
 *
 * TODO: Not Implemented!!!
 *
 * TODO: must finished implementation!!! Special attention
 * TODO: to absorption computation
 *
 * @param layer
 * @param distance
 * @param coh_tmm_data
 * @param position_resolved_data
 * @return
 */
uint8_t position_resolved(
    const uint8_t layer,
    double distance,
    CohTmmData* coh_tmm_data,
    PositionResolvedData* position_resolved_data
)
{
    double complex v;
    double complex w;

    if (layer > 0)
    {
        // TODO: v,w = coh_tmm_data["vw_list"][layer]
        // const double complex v = coh_tmm_data->vw_list[layer][0];
        // const double complex w = coh_tmm_data->vw_list[layer][1];

        // Temporary so I can implement the rest of the function
        v = coh_tmm_data->r;
        w = 1;
    }
    else
    {
        v = 1;
        w = coh_tmm_data->r;
    }
    double complex kz = coh_tmm_data->kz_list[layer];
    const double complex th = coh_tmm_data->th_list[layer];
    const double complex n = coh_tmm_data->n_list[layer];
    const double complex n_0 = coh_tmm_data->n_list[0];
    const double th_0 = coh_tmm_data->th_0;
    const double pol = coh_tmm_data->pol;

    if (
        !(
            layer >= 1 && 0 <= distance <= coh_tmm_data->d_list[layer]
            || layer == 0 && distance <= 0
        )
    )
    {
        return 1;
    }

    // Amplitude of forward-moving wave is Ef, backwards is Eb
    const double complex Ef = v * cexp(1 * I * kz * distance);
    const double complex Eb = w * cexp(-1 * I * kz * distance);

    // Poynting vector
    double complex poyn;
    if (pol == 0)  // s-polarization
    {
        poyn = (
            creal( n * cos(th) * conj(Ef + Eb) * (Ef - Eb) )
            / creal( n_0 * cos(th_0) )
        );
    } else  // p-polarization
    {
        poyn = (
            creal( n * conj( cos(th) ) * (Ef + Eb) * conj( Ef - Eb ) )
            / creal( n_0 * conj( cos(th_0) ) )
        );
    }

    // Absorbed energy density
    double complex absor;
    if (pol == 0)  // s-polarization
    {
        absor = (
            // cimag( n * cos(th) * kz * abs(Ef + Eb) ** 2 )  // make permanent
            cimag( n * cos(th) * kz * abs(Ef + Eb) * 2 )  // temp
            / creal( n_0 * cos(th_0) )
        );
    } else  // p-polarization
    {
        absor = (
            // cimag( n * conj( cos(th) ) * ( kz * abs(Ef - Eb) ** 2 - conj(kz) * abs(Ef + Eb) ** 2 ) )  // make permanent
            cimag( n * conj( cos(th) ) * ( kz * abs(Ef - Eb) * 2 - conj(kz) * abs(Ef + Eb) * 2 ) )  // temp
            / creal( n_0 * conj( cos(th_0) ) )
        );
    }

    // Electric field
    double complex Ex;
    double complex Ey;
    double complex Ez;
    if (pol == 0)  // s-polarization
    {
        Ex = 0.0;
        Ey = Ef + Eb;
        Ez = 0.0;
    } else  // p-polarization
    {
        Ex = (Ef - Eb) * cos(th);
        Ey = 0.0;
        Ez = (-Ef - Eb) * sin(th);
    }

    position_resolved_data->poyn = poyn;
    position_resolved_data->absor = absor;
    position_resolved_data->Ex = Ex;
    position_resolved_data->Ey = Ey;
    position_resolved_data->Ez = Ez;

    return 0;
}


/**
 * @brief
 *
 * Done!!!
 *
 * For large distance, layer = len(d_list), even though d_list[layer] does not
 * exist in this case. For negative distance, return [-1, distance]
 *
 * @param d_list The list of thicknesses of layers, all of which are finite
 * @param d_list_size Size of d_list provided from the outer scope.
 * @param distance
 * @param interface_info Provide array of size two from the outer scope.
 * @return
 */
uint8_t find_in_structure(
    const double d_list[],
    const uint8_t d_list_size,
    double distance,
    double interface_info[]
)
{
    for (int i = 0; i < d_list_size; i++)
    {
        if (d_list[i] == INFINITY)
        {
            printf("this function expects finite arguments\n");
            return 1;  // error
        }
    }
    if (distance < 0.0)
    {
        interface_info[0] = -1.0;
        interface_info[1] = distance;
    } else
    {
        uint8_t layer = 0;
        while ( ( layer < d_list_size ) && ( distance >= d_list[layer] ) )
        {
            distance -= d_list[layer];
            layer++;
        }
        interface_info[0] = layer;
        interface_info[1] = distance;
    }

    return 0;
}


/**
 * @brief
 *
 * Done!!!
 *
 * @param d_list
 * @param d_list_size Size of d_list provided from the outer scope.
 * @param distance
 * @param interface_info Provide array of size two from the outer scope.
 * @return
 */
uint8_t find_in_structure_with_inf(
    const double d_list[],
    const uint8_t d_list_size,
    const double distance,
    double interface_info[]
)
{
    if (distance < 0.0)
    {
        interface_info[0] = 0.0;
        interface_info[1] = distance;
    } else
    {
        find_in_structure(d_list, d_list_size, distance, interface_info);

        const double interface_info_0 = interface_info[0];
        interface_info[0] = interface_info_0 + 1;
    }

    return 0;
}


/**
 * @brief Gives the location of the start of any given layer
 *
 * Done!!!
 *
 * The location of the start of any given layer, relative to the front
 * of the whole multilayer structure. (i.e. the start of layer 1)
 *
 * @param d_list The list of thicknesses of layers
 * @param d_list_size The number of layers
 * @param final_answer
 * @return
 */
uint8_t layer_starts(
    const double d_list[], const uint8_t d_list_size, double final_answer[]
)
{
    // TODO: Confirm that -1 * INFINITY if equivalent to -np.inf
    final_answer[0] = -1 * INFINITY;
    final_answer[1] = 0.0;
    for (int i = 2; i < d_list_size; i++)
    {
        final_answer[i] = final_answer[i - 1] + d_list[i - 1];
    }
    return 0;
}


/**
 * @brief An array listing what proportion of light is absorbed in each layer.
 *
 * Done!!!
 *
 * @param coh_tmm_data
 * @param num_layers
 * @param final_answer
 * @return
 */
uint8_t absorp_in_each_layer(
    CohTmmData* coh_tmm_data,
    const uint8_t num_layers,
    double final_answer[]
)
{
    double power_entering_each_layer[num_layers];
    power_entering_each_layer[0] = 1.0;
    power_entering_each_layer[1] = coh_tmm_data->power_entering;
    power_entering_each_layer[num_layers - 1] = coh_tmm_data->T;

    for (int i = 2; i < num_layers - 1; i++)
    {
        PositionResolvedData position_resolved_data;
        position_resolved(i, 0, coh_tmm_data, &position_resolved_data);
        power_entering_each_layer[i] = position_resolved_data.poyn;
    }

    for (int i = 0; i < num_layers - 1; i++)
    {
        final_answer[i] = (
            -1
            * (power_entering_each_layer[i] - power_entering_each_layer[i + 1])
        );
    }
    final_answer[num_layers - 1] = power_entering_each_layer[num_layers - 1];

    return 0;
}


/**
 * @brief Helper function for inc_tmm. Groups and sorts layer information.
 *
 * See coh_tmm for definitions of n_list, d_list.
 *
 * c_list is "coherency list". Each entry should be 'i' for incoherent or 'c'
 * for 'coherent'.
 *
 * @param n_list
 * @param d_list
 * @param num_layers Number of layers
 * @param c_list coherency list; entries should be 'i' (0) for incoherent or
 *     'c' (1) for 'coherent'
 * @return
 */
uint8_t inc_group_layers(
    double complex n_list[],
    double d_list[],
    uint8_t c_list[],
    const uint8_t num_layers
)
{
    // d_list must start and end with INFINITY
    if (d_list[0] != INFINITY || d_list[num_layers - 1] != INFINITY)
    {
        printf("d_list must start and end with inf!\n");
        // return with error
        return 1;
    }

    // c_list should start and end with 0 (incoherent)
    if (c_list[0] != 0 || c_list[num_layers - 1] != 0)
    {
        printf("c_list should start and end with 0 (incoherent)\n");
        return 1;  // return with error
    }

    uint8_t inc_index = 0;
    uint8_t stack_index = 0;
    bool stack_in_progress = false;

    for (int alllayer_index = 0; alllayer_index < num_layers; alllayer_index++)
    {
        if (c_list[alllayer_index] == 1)  // coherent layer
        {
            printf("coherent layer");
            // TODO: inc_from_all.append(nan)
            if (!stack_in_progress)  // this layer is starting new stack
            {
                stack_in_progress = true;
                uint8_t within_stack_index = 1;
            } else  // another coherent layer in the same stack
            {
                // TODO: within_stack_index += 1
            }

        } else if (c_list[alllayer_index] == 0)  // incoherent layer
        {
            printf("incoherent layer");
            if (!stack_in_progress)  // previous layer was also incoherent
            {
                printf("previous layer was also incoherent");
            } else  // previous layer was coherent
            {
                printf("previous layer was coherent");
            }

        } else
        {
            printf("[ValueError] Error: c_list entries must be 'i' or 'c'!");
            return 1;  // return with error
        }
    }

    return 0;
}


/**
 * @brief Incoherent, or partly-incoherent-partly-coherent, transfer matrix method.
 *
 *
 *
 * @param pol
 * @param n_list
 * @param d_list
 * @param c_list
 * @param num_layers
 * @param th_0
 * @param lam_vac
 * @return
 */
uint8_t inc_tmm(
    uint8_t pol,
    double complex n_list[],
    double d_list[],
    uint8_t c_list[],
    const uint8_t num_layers,
    double th_0,
    double lam_vac
)
{
    // TODO: Input tests
    // inc_group_layers(n_list, d_list, c_list);

    // th_list is a list with, for each layer, the angle that the light
    // travels through the layer. Computed with Snell's law. Note that
    // the "angles" may be complex!
    double complex th_list[num_layers];
    list_snell(n_list, num_layers, th_0, th_list);

    // coh_tmm_data_list[i] is the output of coh_tmm for the i'th stack
    // TODO: coh_tmm_data_list = []
    // coh_tmm_bdata_list[i] is the same stack as coh_tmm_data_list[i] but
    // with order of layers reversed
    // TODO: coh_tmm_bdata_list = []

    // num_stacks placeholder
    uint8_t num_stacks = 4;
    for (int i = 0; i < num_stacks; i++)
    {

    }

    // P_list[i] is fraction not absorbed in a single pass through i'th incoherent
    // layer.
    // num_stacks placeholder
    uint8_t num_inc_layers = 4;
    for (int inc_index = 1; inc_index < num_inc_layers - 1; inc_index++)
    {

    }

    // T_list[i,j] and R_list[i,j] are transmission and reflection powers,
    // respectively, coming from the i'th incoherent layer, going to the j'th
    // incoherent layer. Only need to calculate this when j=i+1 or j=i-1.
    // (2D array is overkill but helps avoid confusion.)
    // initialize these arrays
    double complex t_list[num_inc_layers][num_inc_layers];
    tmm_matrix_zeros(num_inc_layers, t_list);
    double complex r_list[num_inc_layers][num_inc_layers];
    tmm_matrix_zeros(num_inc_layers, r_list);
    for (int inc_index = 0; inc_index < num_inc_layers - 1; inc_index++)  // looking at interface i -> i + 1
    {

    }

    // L is the transfer matrix from the i'th to (i+1)st incoherent layer, see
    // https://arxiv.org/abs/1603.02720
    for (int i = 1; i < num_inc_layers - 1; i++)
    {

    }

    // VW_list[n] = [V_n, W_n], the forward- and backward-moving intensities
    // at the beginning of the n'th incoherent layer. VW_list[0] is undefined
    // because 0'th layer has no beginning.
    for (int i = num_layers - 2; i > 0; i--)
    {

    }

    // stackFB_list[n]=[F,B] means that F is light traveling forward towards n'th
    // stack and B is light traveling backwards towards n'th stack.
    // Reminder: inc_from_stack[i] = j means that the i'th stack comes after the
    // layer with incoherent index j.

    // power_entering_list[i] is the normalized Poynting vector crossing the
    // interface into the i'th incoherent layer from the previous (coherent or
    // incoherent) layer. See https://arxiv.org/abs/1603.02720 .
    for (int i = 1; i < num_inc_layers; i++)
    {

    }

    return 0;
}


/**
 * @brief
 *
 * @param inc_data
 * @return
 */
uint8_t inc_absorp_in_each_layer(double inc_data)
{
    return 0;
}


/**
 * @brief Outputs an absorp_analytic_fn object for a coherent layer
 *     within a partly-incoherent stack.
 *
 * inc_data is output of inc_tmm()
 *
 * @param layer
 * @param inc_data
 * @return
 */
uint8_t inc_find_absorp_analytic_fn(uint8_t layer, double inc_data)
{
    return 0;
}
