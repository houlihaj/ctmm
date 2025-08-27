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
#include "tmm_incoherent.h"
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
        for (int j = 0; j < 2; j++)
        {
            coh_tmm_data->vw_list[i * 2] = vw_list[i][0];  // TODO: needs thorough validation
            coh_tmm_data->vw_list[i * 2 + 1] = vw_list[i][1];  // TODO: needs thorough validation
        }
        coh_tmm_data->kz_list[i] = kz_list[i];
        coh_tmm_data->th_list[i] = th_list[i];
        coh_tmm_data->n_list[i] = n_list[i];
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
 * @brief Calculate the Poynting vector, absorbed energy density,
 * and E-field at a specific location.
 *
 * Done!!!
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
        // v = coh_tmm_data->vw_list[layer][0];
        v = coh_tmm_data->vw_list[layer * 2];  // TODO: needs thorough validation
        // w = coh_tmm_data->vw_list[layer][1];
        w = coh_tmm_data->vw_list[layer * 2 + 1];  // TODO: needs thorough validation
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
            (layer >= 1) && (0 <= distance && distance <= coh_tmm_data->d_list[layer])
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
        double complex Esum = Ef + Eb;
        absor = (
            cimag( n * cos(th) * kz * ( creal(Esum) * creal(Esum) + cimag(Esum) * cimag(Esum) ) )
            / creal( n_0 * cos(th_0) )
        );
    } else  // p-polarization
    {
        double complex Esum = Ef + Eb;
        double complex Ediff = Ef - Eb;
        absor = (
            cimag(
                n * conj(cos(th))
                * (
                    kz
                    * ( creal(Ediff) * creal(Ediff) + cimag(Ediff) * cimag(Ediff) )
                    - conj(kz)
                    * ( creal(Esum) * creal(Esum) + cimag(Esum) * cimag(Esum) )
                )
            )
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
 * TODO: must thoroughly validate implementation!!!
 *
 * See coh_tmm for definitions of n_list, d_list.
 *
 * c_list is "coherency list". Each entry should be 'i' for incoherent or 'c'
 * for 'coherent'.
 *
 * A "stack" is a group of one or more consecutive coherent layers. A "stack
 * index" labels the stacks 0,1,2,.... The "within-stack index" counts the
 * coherent layers within the stack 1,2,3... [index 0 is the incoherent layer
 * before the stack starts]
 *
 * An "incoherent layer index" labels the incoherent layers 0,1,2,...
 *
 * An "alllayer index" labels all layers (all elements of d_list) 0,1,2,...
 *
 * Returns info about how the layers relate:
 *
 * - stack_d_list[i] = list of thicknesses of each coherent layer in the i'th
 *   stack, plus starting and ending with "inf"
 * - stack_n_list[i] = list of refractive index of each coherent layer in the
 *   i'th stack, plus the two surrounding incoherent layers
 * - all_from_inc[i] = j means that the layer with incoherent index i has
 *   alllayer index j
 * - inc_from_all[i] = j means that the layer with alllayer index i has
 *   incoherent index j. If j = nan then the layer is coherent.
 * - all_from_stack[i1][i2] = j means that the layer with stack index i1 and
 *   within-stack index i2 has alllayer index j
 * - stack_from_all[i] = [j1 j2] means that the layer with alllayer index i is
 *   part of stack j1 with withinstack-index j2. If stack_from_all[i] = nan
 *   then the layer is incoherent
 * - inc_from_stack[i] = j means that the i'th stack comes after the layer
 *   with incoherent index j, and before the layer with incoherent index j+1.
 * - stack_from_inc[i] = j means that the layer with incoherent index i comes
 *   immediately after the j'th stack. If j=nan, it is not immediately
 *   following a stack.
 *
 * @param inc_group_layers_data
 * @param n_list
 * @param d_list
 * @param c_list coherency list; entries should be 'i' (0) for incoherent or
 *     'c' (1) for 'coherent'
 * @param num_layers Number of layers; same as size of n_list and d_list
 * @param num_inc_layers Number of incoherent layers
 * @param num_stacks Number of stacks
 * @return
 */
uint8_t inc_group_layers(
    IncGroupLayersData* inc_group_layers_data,
    double complex n_list[],
    double d_list[],
    uint8_t c_list[],
    const uint8_t num_layers,
    const uint8_t num_inc_layers,
    const uint8_t num_stacks
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

    // // count the number of incoherent layers in c_list
    // uint8_t num_inc_layers = 0;  // number of incoherent layers
    // for (int i = 0; i < num_layers; i++)
    // {
    //     if (c_list[i] == 0)  // incoherent layer
    //     {
    //         num_inc_layers++;
    //     }
    // }
    //
    // // count the number of stacks in c_list
    // uint8_t num_stacks = 0;  // number of stacks
    // for (int i = 1; i < num_layers; i++)
    // {
    //     if (c_list[i - 1] == 0 && c_list[i] == 1)
    //     {
    //         num_stacks++;
    //     }
    // }

    uint8_t inc_index = 0;
    uint8_t stack_index = 0;

    double* stack_d_list = inc_group_layers_data->stack_d_list;
    double complex* stack_n_list = inc_group_layers_data->stack_n_list;
    uint8_t* all_from_inc = inc_group_layers_data->all_from_inc;
    uint8_t* inc_from_all = inc_group_layers_data->inc_from_all;
    uint8_t* all_from_stack = inc_group_layers_data->all_from_stack;
    uint8_t* stack_from_all = inc_group_layers_data->stack_from_all;
    uint8_t* inc_from_stack = inc_group_layers_data->inc_from_stack;
    uint8_t* stack_from_inc = inc_group_layers_data->stack_from_inc;

    const uint8_t nan = 255;  // use in place of np.nan
    uint8_t within_stack_index = 0;  // TODO: must validate!!! differs from TMM

    double ongoing_stack_d_list[num_stacks * 2];
    double complex ongoing_stack_n_list[num_stacks * 2];

    bool stack_in_progress = false;

    for (int alllayer_index = 0; alllayer_index < num_layers; alllayer_index++)
    {
        if (c_list[alllayer_index] == 1)  // coherent layer
        {
            printf("coherent layer\n");
            inc_from_all[alllayer_index] = nan;  // gtg; array size is equivalent to num_layers
            if (!stack_in_progress)  // this layer is starting new stack
            {
                within_stack_index = 0;
                stack_in_progress = true;
                ongoing_stack_d_list[within_stack_index] = INFINITY;  // gtg
                ongoing_stack_d_list[within_stack_index + 1] = d_list[alllayer_index];  // gtg
                ongoing_stack_n_list[within_stack_index] = n_list[alllayer_index - 1];  // gtg
                ongoing_stack_n_list[within_stack_index + 1] = n_list[alllayer_index];  // gtg
                stack_from_all[alllayer_index * 2] = stack_index;  // gtg; size is equivalent to (2 * num_coh_layers + 1 * num_inc_layers)
                stack_from_all[alllayer_index * 2 + 1] = 1;  // gtg; size is equivalent to (2 * num_coh_layers + 1 * num_inc_layers)
                all_from_stack[within_stack_index] = alllayer_index - 1;  // gtg
                all_from_stack[within_stack_index + 1] = alllayer_index;  // gtg
                inc_from_stack[stack_index] = inc_index - 1;  // gtg; size is equivalent to num_stacks
                within_stack_index = 1;  // or within_stack_index++
            } else  // another coherent layer in the same stack
            {
                ongoing_stack_d_list[within_stack_index] = d_list[alllayer_index];  // gtg
                ongoing_stack_n_list[within_stack_index] = n_list[alllayer_index];  // gtg
                within_stack_index++;
                stack_from_all[alllayer_index * 2] = stack_index;  // gtg; size is equivalent to (2 * num_coh_layers + 1 * num_inc_layers)
                stack_from_all[alllayer_index * 2 + 1] = within_stack_index;  // gtg; size is equivalent to (2 * num_coh_layers + 1 * num_inc_layers)
                all_from_stack[within_stack_index - 1] = alllayer_index;  // gtg
            }

        } else if (c_list[alllayer_index] == 0)  // incoherent layer
        {
            printf("incoherent layer\n");
            stack_from_all[alllayer_index] = nan;  // gtg; size is equivalent to (2 * num_coh_layers + 1 * num_inc_layers)
            inc_from_all[alllayer_index] = inc_index;  // gtg; array size is equivalent to num_layers
            all_from_inc[inc_index] = alllayer_index;  // gtg; array size is equivalent to num_inc_layers
            if (!stack_in_progress)  // previous layer was also incoherent
            {
                printf("previous layer was also incoherent\n");
                stack_from_inc[inc_index] = nan;  // gtg; size is equivalent to num_inc_layers
            } else  // previous layer was coherent
            {
                printf("previous layer was coherent\n");
                stack_in_progress = false;
                stack_from_inc[inc_index] = stack_index;  // gtg; size is equivalent to num_inc_layers
                ongoing_stack_d_list[within_stack_index] = INFINITY;  // sketchy but gtg
                for (int i = 0; i < stack_index * 2; i++)  // i should not start at zero always
                {
                    stack_d_list[i] = ongoing_stack_d_list[i];
                }
                ongoing_stack_n_list[within_stack_index] = n_list[alllayer_index];  // sketchy but gtg
                for (int i = 0; i < stack_index * 2; i++)  // i should not start at zero always
                {
                    stack_n_list[i] = ongoing_stack_n_list[i];
                }
                all_from_stack[within_stack_index - 1] = alllayer_index;  // gtg
                stack_index++;
            }

            inc_index++;

        } else
        {
            printf("[ValueError] Error: c_list entries must be 'i' or 'c'!\n");
            return 1;  // return with error
        }
    }

    inc_group_layers_data->stack_d_list = stack_d_list;
    inc_group_layers_data->stack_n_list = stack_n_list;
    inc_group_layers_data->all_from_inc = all_from_inc;
    inc_group_layers_data->inc_from_all = inc_from_all;
    inc_group_layers_data->all_from_stack = all_from_stack;
    inc_group_layers_data->stack_from_all = stack_from_all;
    inc_group_layers_data->inc_from_stack = inc_from_stack;
    inc_group_layers_data->stack_from_inc = stack_from_inc;
    inc_group_layers_data->num_stacks = num_stacks;
    inc_group_layers_data->num_inc_layers = num_inc_layers;
    inc_group_layers_data->num_layers = num_layers;

    return 0;
}


/**
 * @brief Incoherent, or partly-incoherent-partly-coherent, transfer matrix method.
 *
 * TODO: must implement!!!
 *
 * @param pol
 * @param n_list
 * @param d_list
 * @param c_list
 * @param num_layers
 * @param num_inc_layers
 * @param num_stacks
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
    const uint8_t num_inc_layers,
    const uint8_t num_stacks,
    double th_0,
    double lam_vac
)
{
    // TODO: Input tests

    // Get the incoherent group layers data
    IncGroupLayersData inc_group_layers_data;
    IncGroupLayersData_create(
        &inc_group_layers_data, num_layers, num_inc_layers, num_stacks
    );
    inc_group_layers(
        &inc_group_layers_data, n_list, d_list, c_list, num_layers, num_inc_layers, num_stacks
    );
    double complex* stack_n_list = inc_group_layers_data.stack_n_list;  // TODO: validate member access is correct
    double* stack_d_list = inc_group_layers_data.stack_d_list;
    uint8_t* all_from_stack = inc_group_layers_data.all_from_stack;
    uint8_t* all_from_inc = inc_group_layers_data.all_from_inc;
    uint8_t* stack_from_inc = inc_group_layers_data.stack_from_inc;
    uint8_t* inc_from_stack = inc_group_layers_data.inc_from_stack;

    // th_list is a list with, for each layer, the angle that the light
    // travels through the layer. Computed with Snell's law. Note that
    // the "angles" may be complex!
    double complex th_list[num_layers];
    list_snell(n_list, num_layers, th_0, th_list);

    // coh_tmm_data_list[i] is the output of coh_tmm for the i'th stack
    CohTmmData coh_tmm_data_list[num_stacks];
    // coh_tmm_bdata_list[i] is the same stack as coh_tmm_data_list[i] but
    // with order of layers reversed
    CohTmmData coh_tmm_bdata_list[num_stacks];

    for (int i = 0; i < num_stacks; i++)
    {
        // CohTmmData coh_tmm_data;  // allocate struct in outer scope
        // const uint8_t num_coh_layers = (num_layers - num_inc_layers);
        // CohTmmData_create(&coh_tmm_data, num_layers);
        // coh_tmm(
        //     pol,  // const uint8_t pol,
        //     stack_n_list[i],  // const double complex n_list[],
        //     stack_d_list[i],  // const double d_list[],
        //     num_coh_layers,  // const uint8_t num_layers,
        //     th_list[all_from_stack[i][0]],  // const double th_0,
        //     lam_vac,  // const double lam_vac,
        //     &coh_tmm_data  // CohTmmData* coh_tmm_data
        // );
        //
        // coh_tmm_data_list[i] = coh_tmm_data;
    }

    // P_list[i] is fraction not absorbed in a single pass through i'th incoherent
    // layer.
    double P_list[num_inc_layers];
    for (int inc_index = 1; inc_index < num_inc_layers - 1; inc_index++)  // skip 0'th and last (infinite)
    {
        uint8_t i = all_from_inc[inc_index];
        P_list[inc_index] = (
            exp( -4 * PI * d_list[i] * cimag(n_list[i] * cos(th_list[i])) / lam_vac )
        );
        // For a very opaque layer, reset P to avoid divide-by-0 and similar
        // errors.
        if (P_list[inc_index] < 1.0e-30)
        {
            P_list[inc_index] = 1.0e-30;
        }
    }

    // T_list[i,j] and R_list[i,j] are transmission and reflection powers,
    // respectively, coming from the i'th incoherent layer, going to the j'th
    // incoherent layer. Only need to calculate this when j=i+1 or j=i-1.
    // (2D array is overkill but helps avoid confusion.)
    // initialize these arrays
    double complex T_list[num_inc_layers][num_inc_layers];
    tmm_matrix_zeros(num_inc_layers, T_list);
    double complex R_list[num_inc_layers][num_inc_layers];
    tmm_matrix_zeros(num_inc_layers, R_list);
    for (int inc_index = 0; inc_index < num_inc_layers - 1; inc_index++)  // looking at interface i -> i + 1
    {
        uint8_t alllayer_index = all_from_inc[inc_index];
        uint8_t nextstack_index = stack_from_inc[inc_index + 1];
        const uint8_t nan = 255;  // use in place of np.nan
        if (nextstack_index == nan)  // next layer is incoherent
        {
            double R;
            interface_R(
                pol,  // const uint8_t polarization,
                n_list[alllayer_index],  // const double complex n_i,
                n_list[alllayer_index + 1],  // const double complex n_f,
                th_list[alllayer_index],  // const double th_i,
                th_list[alllayer_index + 1],  // const double th_f,
                &R  // double* R
            );
            R_list[inc_index][inc_index + 1] = R;

            double T;
            interface_T(
                pol,  // const uint8_t polarization,
                n_list[alllayer_index],  // const double complex n_i,
                n_list[alllayer_index + 1],  // const double complex n_f,
                th_list[alllayer_index],  // const double th_i,
                th_list[alllayer_index + 1],  // const double th_f,
                &T  // double* T
            );
            T_list[inc_index][inc_index + 1] = T;

            interface_R(
                pol,  // const uint8_t polarization,
                n_list[alllayer_index + 1],  // const double complex n_i,
                n_list[alllayer_index],  // const double complex n_f,
                th_list[alllayer_index + 1],  // const double th_i,
                th_list[alllayer_index],  // const double th_f,
                &R  // double* R
            );
            R_list[inc_index + 1][inc_index] = R;

            interface_T(
                pol,  // const uint8_t polarization,
                n_list[alllayer_index + 1],  // const double complex n_i,
                n_list[alllayer_index],  // const double complex n_f,
                th_list[alllayer_index + 1],  // const double th_i,
                th_list[alllayer_index],  // const double th_f,
                &T  // double* T
            );
            T_list[inc_index + 1][inc_index] = T;
        } else  // next layer is coherent
        {
            R_list[inc_index][inc_index + 1] = (
                coh_tmm_data_list[nextstack_index].R
            );
            T_list[inc_index][inc_index + 1] = (
                coh_tmm_data_list[nextstack_index].T
            );

            R_list[inc_index + 1][inc_index] = (
                coh_tmm_bdata_list[nextstack_index].R
            );
            T_list[inc_index + 1][inc_index] = (
                coh_tmm_bdata_list[nextstack_index].T
            );
        }
    }

    // L is the transfer matrix from the i'th to (i+1)st incoherent layer, see
    // https://arxiv.org/abs/1603.02720
    // TODO: L_list = [nan] # L_0 is not defined because 0'th layer has no beginning.
    // TODO: Ltilde declaration and assignment
    for (int i = 1; i < num_inc_layers - 1; i++)
    {
        double mat1[2][2];
        double mat2[2][2];
        // tmm_scalar_division(mat2, T_list[i][i + 1]);
        // tmm_matrix_product(mat1, mat2);
    }
    double T = 0.0;
    double R = 0.0;

    // VW_list[n] = [V_n, W_n], the forward- and backward-moving intensities
    // at the beginning of the n'th incoherent layer. VW_list[0] is undefined
    // because 0'th layer has no beginning.
    double VW_list[num_inc_layers][2];
    const uint8_t nan = 255;  // use in place of np.nan
    VW_list[0][0] = nan;  // TODO: validate proper value for nan
    VW_list[0][1] = nan;  // TODO: validate proper value for nan
    // 2x1 array (2 rows, 1 column)
    double VW[2][1] = {{T},{0.0}};

    // // 1x2 array (1 rows, 2 column)
    // double complex VW_tr[1][2];
    // tmm_transpose(VW, VW_tr);
    // VW_list[num_layers - 1][0] = VW_tr[0][0];
    // VW_list[num_layers - 1][1] = VW_tr[0][1];

    // for (int i = num_layers - 2; i > 0; i--)  // temporary to snap to previous implementation
    for (int i = num_inc_layers - 2; i > 0; i--)
    {
        // TODO: VW = np.dot(L_list[i], VW)
        // TODO: VW_list[i,:] = np.transpose(VW)

        // tmm_matrix_by_vector(L_list[i], VW, VW);
        // tmm_transpose(VW, VW_tr);
        // VW_list[i][0] = VW_tr[0][0];
        // VW_list[i][1] = VW_tr[0][1];
    }

    // stackFB_list[n]=[F,B] means that F is light traveling forward towards n'th
    // stack and B is light traveling backwards towards n'th stack.
    // Reminder: inc_from_stack[i] = j means that the i'th stack comes after the
    // layer with incoherent index j.
    // TODO: stackFB_list = []
    double complex* stackFB_list = inc_group_layers_data.inc_tmm_data.stackFB_list;
    for (int stack_i = 0; stack_i < num_inc_layers; stack_i++)  // TODO: must confirm inc_from_stack length in tmm_core.py
    {
        double F;
        uint8_t prev_inc_index = inc_from_stack[stack_i];  // TODO: must validate!!!
        if (prev_inc_index == 0)  // stack starts right after semi-infinite layer.
        {
            F = 1.0;
        } else
        {
            F = VW_list[prev_inc_index][0] * P_list[prev_inc_index];
        }
        double B = VW_list[prev_inc_index + 1][1];
        // TODO: stackFB_list.append([F,B])
        stackFB_list[stack_i * 2] = F;  // TODO: must validate!!!
        stackFB_list[stack_i * 2 + 1] = B;  // TODO: must validate!!!
    }

    // power_entering_list[i] is the normalized Poynting vector crossing the
    // interface into the i'th incoherent layer from the previous (coherent or
    // incoherent) layer. See https://arxiv.org/abs/1603.02720 .
    double power_entering_list[num_inc_layers];  // TODO: update with struct member
    power_entering_list[0] = 1.0;  // "1" by convention for infinite 0th layer.
    for (int i = 1; i < num_inc_layers; i++)
    {
        uint8_t prev_stack_index = stack_from_inc[1];
    }

    IncGroupLayersData_destroy(&inc_group_layers_data);

    return 0;
}


/**
 * @brief
 *
 * TODO: must implement!!!
 *
 * @param inc_tmm_data
 * @param absorp_list
 * @return
 */
uint8_t inc_absorp_in_each_layer(
    // IncTmmData* inc_tmm_data, double absorp_list[]
    // TODO: confirm this is the proper approach
    IncGroupLayersData* inc_group_layers_data,
    IncTmmData* inc_tmm_data,
    double absorp_list[]
)
{
    // Reminder: inc_from_stack[i] = j means that the i'th stack comes after the
    // layer with incoherent index j.
    // Reminder: stack_from_inc[i] = j means that the layer
    // with incoherent index i comes immediately after the j'th stack (or j=nan
    // if it's not immediately following a stack).

    const uint8_t num_inc_layers = inc_tmm_data->num_inc_layers;

    uint8_t* stack_from_inc = inc_group_layers_data->stack_from_inc;  // TODO: must add member to struct
    double* power_entering_list = inc_tmm_data->power_entering_list;
    // stackFB_list[n]=[F,B] means that F is light traveling forward towards n'th
    // stack and B is light traveling backwards towards n'th stack.
    double complex* stackFB_list = inc_tmm_data->stackFB_list;

    for (int i = 0; i < num_inc_layers - 1; i++)
    {
        const uint8_t nan = 255;  // use in place of np.nan
        if (stack_from_inc[i + 1] == nan)
        {
            // case that incoherent layer i is right before another incoherent layer
            absorp_list[i] = power_entering_list[i] - power_entering_list[i + 1];
        } else  // incoherent layer i is immediately before a coherent stack
        {
            uint8_t j = stack_from_inc[i + 1];
            CohTmmData coh_tmm_data = inc_tmm_data->coh_tmm_data_list[j];
            CohTmmData coh_tmm_bdata = inc_tmm_data->coh_tmm_bdata_list[j];
            // First, power in the incoherent layer...
            // TODO: temporary!!!
            absorp_list[i] = power_entering_list[i] - power_entering_list[i + 1];
            // Next, power in the coherent stack...
        }
    }

    // final semi-infinite layer
    absorp_list[num_inc_layers - 1] = inc_tmm_data->T;

    return 0;
}


/**
 * @brief Outputs an absorp_analytic_fn object for a coherent layer
 *     within a partly-incoherent stack.
 *
 * inc_data is output of inc_tmm()
 *
 * @param layer
 * @param inc_tmm_data
 * @return
 */
uint8_t inc_find_absorp_analytic_fn(
    uint8_t layer,
    // IncTmmData* inc_tmm_data
    IncGroupLayersData* inc_group_layers_data
)
{
    // uint8_t* j = inc_group_layers_data->stack_from_all[layer];
    uint8_t* j = inc_group_layers_data->stack_from_all;

    // TODO: add is coherent layer check; refer to tmm_core.py
    // layer must be coherent for this function!

    uint8_t stackindex = j[layer * 2];  // TODO: must validate!!!
    uint8_t withinstackindex = j[layer * 2 + 1];  // TODO: must validate!!!

    AbsorpAnalyticFn forward_absorp_fn;
    // CohTmmData coh_tmm_data;
    CohTmmData coh_tmm_data = (
        inc_group_layers_data->inc_tmm_data.coh_tmm_data_list[stackindex]
    );
    fill_in(
        &forward_absorp_fn, &coh_tmm_data, withinstackindex
    );
    // scale(&forward_absorp_fn, );
    flip(&forward_absorp_fn);

    AbsorpAnalyticFn back_absorp_fn;
    CohTmmData coh_tmm_bdata = (
        inc_group_layers_data->inc_tmm_data.coh_tmm_bdata_list[stackindex]
    );
    fill_in(&back_absorp_fn, &coh_tmm_bdata, layer);
    // scale(&back_absorp_fn, );
    flip(&back_absorp_fn);

    add(&forward_absorp_fn, &back_absorp_fn);

    return 0;
}
