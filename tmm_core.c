//
// Created by johnh on 8/2/2025.
//

#include <stdint.h>
#include <complex.h>
#include <math.h>
#include "tmm_core.h"
#include "tmm_absorp_fn.h"


/**
 * @brief asdfas
 *
 * This struct holds integer values for the x and y coordinates.
 */
typedef struct {
    double r;  // complex reflectivity amplitude
    double t;  // complex transmissivity amplitude
    double R;
    double T;
    double power_entering;
    double vw_list;
    double kz_list;
    double th_list;
    double pol;
    double n_list;
    double d_list;
    double th_0;
    double lam_vac;
} CohTmmData;


typedef struct {
    double psi;  // units of radians
    double delta;  // units of radians
} EllipsData;


/**
 * @brief Makes a 2x2 array of [[a,b],[c,d]]
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
uint8_t make_2x2_array(
    const double a,
    const double b,
    const double c,
    const double d,
    double complex matrix[2][2]
)
{
    matrix[0][0] = a;
    matrix[0][1] = b;
    matrix[1][0] = c;
    matrix[1][1] = d;

    return 0;
}


/**
 * @brief Compute whether-or-not this is the forward-traveling wave
 *
 * @param n Medium with complex refractive index, n
 * @param theta Angle of incidence
 * @return
 */
uint8_t is_forward_angle(double complex n, double theta)
{
    double complex ncostheta = n * cos(theta);

    return 0;
}


/**
 * @brief Return the angle theta in layer 2 with refractive index n_2
 *
 * @param n_1
 * @param n_2
 * @param th_1
 * @return
 */
double snell(const double n_1, const double n_2, const double th_1)
{
    const double th_2_guess = asin(n_1 * sin(th_1) / n_2);

    // is_forward_angle()

    return th_2_guess;
}


/**
 * @brief
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
    double angles[]
)
{
    const double n_list_0 = n_list[0];
    for (uint8_t i = 0; i < n_list_size; i++)
    {
        angles[i] = asin( n_list_0 * sin(th_0) / n_list[i] );
    }
    return 0;
}


/**
 * @brief Complex reflection amplitude from Fresnel equations
 *
 * @param polarization
 * @param n_i
 * @param n_f
 * @param th_i
 * @param th_f
 * @return
 */
double complex interface_r(
    const uint8_t polarization, const double complex n_i, const double complex n_f, const double th_i, const double th_f
)
{
    double complex reflectivity;

    if (polarization == 0)  // s-polarization
    {
        reflectivity = (
            ( n_i * cos(th_i) - n_f * cos(th_f) )
            / ( n_i * cos(th_i) + n_f * cos(th_f) )
        );
    } else  // p-polarization
    {
        reflectivity = (
            ( n_f * cos(th_i) - n_i * cos(th_f) )
            / ( n_f * cos(th_i) + n_i * cos(th_f) )
        );
    }

    return reflectivity;
}


/**
 * @brief
 *
 *
 *
 * @param
 * @param
 * @return
 */
double complex interface_t(
    uint8_t polarization, double complex n_i, double complex n_f, double th_i, double th_f
)
{
    double complex transmissivity;

    if (polarization == 0)  // s-polarization
    {
        transmissivity = (
            2 * n_i * cos(th_i)
            / ( n_i * cos(th_i) + n_f * cos(th_f) )
        );
    } else  // p-polarization
    {
        transmissivity = (
            2 * n_i * cos(th_i)
            / ( n_f * cos(th_i) + n_i * cos(th_f) )
        );
    }

    return transmissivity;
}


/**
 * @brief Calculate reflected power R, starting with reflection amplitude r.
 *
 *
 *
 * @param
 * @param
 * @return
 */
double R_from_r(const double complex r)
{
    // Method 1: Manually compute the modulus squared
    const double real_part = creal(r);
    const double imag_part = cimag(r);
    const double reflectance = real_part * real_part + imag_part * imag_part;

    return reflectance;
}


/**
 * @brief Calculate transmitted power T, starting with transmission amplitude t.
 *
 *
 *
 * @param pol
 * @param t
 * @param n_i
 * @param n_f
 * @param th_i
 * @param th_f
 * @return
 */
double T_from_t(
    const uint8_t pol,
    const double complex t,
    const double complex n_i,
    const double complex n_f,
    const double th_i,
    const double th_f
)
{
    // Method 1: Manually compute the modulus squared
    const double real_part = creal(t);
    const double imag_part = cimag(t);

    double transmittance;

    if (pol == 0)  // s-polarization
    {
        transmittance = (
            (real_part * real_part + imag_part * imag_part)
            * ( creal( n_f * cos(th_f) ) / creal( n_i * cos(th_i) ) )
        );
    } else  // p-polarization
    {
        transmittance = (
            (real_part * real_part + imag_part * imag_part)
            * ( creal( n_f * conj(cos(th_f)) ) / creal( n_i * conj(cos(th_i)) ) )
        );
    }

    return transmittance;
}


/**
 * @brief Calculate the power entering the first interface, starting with
 * complex reflection amplitude, r.
 *
 *
 *
 * @param pol
 * @param r
 * @param n_i The refractive index of incident medium
 * @param th_i The complex propogation angle through incident medium
 * @return
 */
double power_entering_from_r(
    const uint8_t pol,
    const double complex r,
    const double complex n_i,
    const double th_i
)
{
    double power;

    if (pol == 0)  // s-polarization
    {
        power = (
            creal( n_i * cos(th_i) * (1 + conj(r)) * (1 - r) )
            / creal( n_i * cos(th_i) )
        );
    } else  // p-polarization
    {
        power = (
            creal( n_i * conj(cos(th_i)) * (1 + r) * (1 - conj(r)) )
            / creal( n_i * conj(cos(th_i)) )
        );
    }

    return power;
}


/**
 * @brief
 *
 *
 *
 * @param
 * @param
 * @return
 */
double interface_R(
    const uint8_t polarization,
    const double complex n_i,
    const double complex n_f,
    const double th_i,
    const double th_f
)
{
    const double complex r = (
        interface_r(polarization, n_i, n_f, th_i, th_f)
    );
    const double R = R_from_r(r);
    return R;
}


/**
 * @brief
 *
 * @param
 * @param
 * @return
 */
double interface_T(
    const uint8_t polarization,
    const double complex n_i,
    const double complex n_f,
    const double th_i,
    const double th_f
)
{
    const double complex t = (
        interface_t(polarization, n_i, n_f, th_i, th_f)
    );
    const double T = T_from_t(polarization, t, n_i, n_f, th_i, th_f);
    return T;
}


/**
 * @brief Main "coherent transfer matrix method" computation.
 *
 *
 *
 * @param pol Light polarization, "s" or "p"
 * @param n_list
 * @param d_list
 * @param num_layers
 * @param th_0 The angle of incidence
 * @param lam_vac The vacuum wavelength of the light
 * @return
 */
uint8_t coh_tmm(
    const uint8_t pol,
    const double complex n_list[],
    const double d_list[],
    const uint8_t num_layers,
    const double th_0,
    const double lam_vac
)
{
    // is_forward_angle();

    double th_list[num_layers];
    list_snell(n_list, num_layers, th_0, th_list);

    // kz is the z-component of (complex) angular wavevector for forward-moving
    // wave. Positive imaginary part means decaying.
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
        delta[i] = kz_list[i] * d_list[i];
    }

    // For a very opaque layer, reset delta to avoid divide-by-0 and similar
    // errors. The criterion imag(delta) > 35 corresponds to single-pass
    // transmission < 1e-30 --- small enough that the exact value doesn't
    // matter.
    for (int i = 1; i < num_layers - 1; i++)
    {
        if ( cimag(delta[i]) > 35 )
        {
            delta[i] = creal(delta[i]) + cimag(35);
        }

        // TODO: Must add opacity check and warning!!!
    }

    // t_list[i,j] and r_list[i,j] are transmission and reflection amplitudes,
    // respectively, coming from i, going to j. Only need to calculate this when
    // j=i+1. (2D array is overkill but helps avoid confusion.)
    double complex t_list[num_layers][num_layers];
    double complex r_list[num_layers][num_layers];
    for (int i = 0; i < num_layers - 1; i++)
    {
        t_list[i][i + 1] = (
            interface_t(
                pol, n_list[i], n_list[i + 1], th_list[i], th_list[i + 1]
            )
        );

        r_list[i][i + 1] = (
            interface_r(
                pol, n_list[i], n_list[i + 1], th_list[i], th_list[i + 1]
            )
        );
    }

    // At the interface between the (n-1)st and nth material, let v_n be the
    // amplitude of the wave on the nth side heading forwards (away from the
    // boundary), and let w_n be the amplitude on the nth side heading backwards
    // (towards the boundary). Then (v_n,w_n) = M_n (v_{n+1},w_{n+1}). M_n is
    // M_list[n]. M_0 and M_{num_layers-1} are not defined.
    // My M is a bit different than Sernelius's, but Mtilde is the same.
    double complex M_list[num_layers][2][2];
    for (int i = 1; i < num_layers - 1; i++)
    {
        continue;
    }
    double complex Mtilde[2][2];
    make_2x2_array(1, 0, 0, 1, Mtilde);

    // Net complex transmission and reflection amplitudes
    double complex r = Mtilde[1][0] / Mtilde[0][0];
    double complex t = 1 / Mtilde[0][0];

    // vw_list[n] = [v_n, w_n]. v_0 and w_0 are undefined because the 0th medium
    // has no left interface.
    double complex vw_list[num_layers][2];
    // TODO: vw = array([[t],[0]])
    // TODO: vw_list[-1,:] = np.transpose(vw)
    // TODO: must confirm loop is implemented correctly
    for (int i = num_layers - 2; i > 0; i--)
    {
        // TODO: vw = np.dot(M_list[i], vw)
        // TODO: vw_list[i,:] = np.transpose(vw)
        continue;
    }

    // Net transmitted and reflected power, as a proportion of the
    // incoming light power.
    const double R = R_from_r(r);
    // TODO: confirm n_list[-1] and n_list[num_layers] are equivalent
    // TODO: confirm th_list[-1] and th_list[num_layers] are equivalent
    const double T = (
        T_from_t(
            pol, t, n_list[0], n_list[num_layers], th_list[0], th_list[num_layers]
        )
    );
    const double power_entering = power_entering_from_r(pol, r, n_list[0], th_0);

    return 0;
}


/**
 * @brief Reverses the order of the stack and then runs coh_tmm()
 *
 *
 *
 * @param
 * @param
 * @return
 */
double coh_tmm_reverse(
    const uint8_t pol,
    const double n_list[],
    const double d_list[],
    const uint8_t num_layers,
    const double th_0,
    const double lam_vac
)
{
    // TODO: confirm n_list[-1] and n_list[num_layers] is equivalent
    const double th_f = snell(n_list[0], n_list[num_layers], th_0);

    double complex n_list_rev[num_layers];
    double d_list_rev[num_layers];
    for (int i = 0; i < num_layers; i++)
    {
        n_list_rev[i] = n_list[num_layers - i];
        d_list_rev[i] = d_list[num_layers - i];
    }

    coh_tmm(pol, n_list_rev, d_list_rev, num_layers, th_f, lam_vac);

    return 0.0;
}


/**
 * @brief Calculates ellipsometric parameters, in radians.
 *
 * Warning: Conventions differ. You may need to subtract pi/2 or whatever.
 *
 * @param n_list
 * @param d_list
 * @param num_layers
 * @param th_0
 * @param lam_vac
 * @return
 */
double ellips(
    const double n_list[],
    const double d_list[],
    const uint8_t num_layers,
    const double th_0,
    const double lam_vac
)
{
    uint8_t s_data = coh_tmm(0, n_list, d_list, num_layers, th_0, lam_vac);
    uint8_t p_data = coh_tmm(1, n_list, d_list, num_layers, th_0, lam_vac);
    // TODO: finish implementation!!!
    return 0.0;
}


/**
 * @brief Calculates reflected and transmitted power for unpolarized light.
 *
 * @param n_list
 * @param d_list
 * @param num_layers
 * @param th_0
 * @param lam_vac
 * @return
 */
double unpolarized_RT(
    const double n_list[],
    const double d_list[],
    const uint8_t num_layers,
    const double th_0,
    const double lam_vac
)
{
    // TODO: finish implementation!!!
    return 0.0;
}


/**
 * @brief Calculate the Poynting vector, absorbed energy density, and E-field at a specific location.
 *
 * @param layer
 * @param distance
 * @param coh_tmm_data
 * @return
 */
double position_resolved(
    uint8_t layer, double distance, double coh_tmm_data
)
{
    return 0.0;
}


/**
 * @brief
 *
 * @param d_list The list of thicknesses of layers, all of which are finite
 * @param d_list_size
 * @param distance
 * @return
 */
uint8_t find_in_structure(
    const double d_list[], const uint8_t d_list_size, const double distance
)
{
    return 0;
}


/**
 * @brief
 *
 * @param d_list
 * @param distance
 * @return
 */
uint8_t find_in_structure_with_inf(
    const double d_list[], const double distance
)
{
    // find_in_structure(d_list, d_list_size, distance);
    return 0;
}


/**
 * @brief Gives the location of the start of any given layer
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
    final_answer[0] = 0.0;  // TODO: replace with -inf
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
 * TODO: must finished implementation!!!
 *
 * @param coh_tmm_data
 * @param num_layers
 * @param final_answer
 * @return
 */
double absorp_in_each_layer(
    double coh_tmm_data,
    const uint8_t num_layers,
    double final_answer[]
)
{
    double power_entering_each_layer[num_layers];
    power_entering_each_layer[0] = 1.0;

    for (int i = 2; i < num_layers - 1; i++)
    {
        // TODO: must finished position_resolved() implementation
        // power_entering_each_layer[i] = position_resolved(i, 0, coh_tmm_data);
        continue;
    }

    // double final_answer[num_layers];
    for (int i = 0; i < num_layers - 1; i++)
    {
        final_answer[i] = power_entering_each_layer[i] - power_entering_each_layer[i + 1];
    }
    final_answer[num_layers] = power_entering_each_layer[num_layers];

    return 0.0;
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
 * @param c_list
 * @return
 */
double inc_group_layers(double n_list[], double d_list[], double c_list[])
{
    return 0.0;
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
 * @param th_0
 * @param lam_vac
 * @return
 */
double inc_tmm(
    uint8_t pol, double n_list[], double d_list[], double c_list[], double th_0, double lam_vac
)
{
    return 0.0;
}


/**
 * @brief
 *
 * @param inc_data
 * @return
 */
double inc_absorp_in_each_layer(double inc_data)
{
    return 0.0;
}


/**
 * @brief Outputs an absorp_analytic_fn object for a coherent layer within a
 *     partly-incoherent stack.
 *
 * inc_data is output of inc_tmm()
 *
 * @param layer
 * @param inc_data
 * @return
 */
double inc_find_absorp_analytic_fn(uint8_t layer, double inc_data)
{
    return 0.0;
}
