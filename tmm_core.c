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
 * @param
 * @param
 * @return
 */
uint8_t make_2x2_array(
    const double a, const double b, const double c, const double d
)
{
    // Must declare the matrix in the outer scope and pass it into
    // this function for initialization of the elements

    // Declare and initialize a 2x2 matrix (example to be used in outer scope)
    double matrix[2][2] = {
        {a, b},
        {c, d}
    };

    // The core of the function is below

    // matrix[0][0] = a
    // matrix[0][1] = b
    // matrix[1][0] = c
    // matrix[1][1] = d

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
    const double n_list[], const uint8_t n_list_size, const double th_0, double angles[]
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
double interface_r(
    const uint8_t polarization, const double n_i, const double n_f, const double th_i, const double th_f
)
{
    double reflectivity;

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
double interface_t(
    uint8_t polarization, double n_i, double n_f, double th_i, double th_f
)
{
    return 0.0;
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
double R_from_r(const double r)
{
    // Must take modulus squared of r
    const double reflectance = r * r;  // temporary

    // // Method 1: Manually compute the modulus squared
    // double real_part = creal(z);
    // double imag_part = cimag(z);
    // double mod_squared = real_part * real_part + imag_part * imag_part;

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
    const uint8_t pol, const double t, const double n_i, const double n_f, const double th_i, const double th_f
)
{
    // Must take modulus squared of t
    const double transmittance = t * t;  // temporary

    // // Method 1: Manually compute the modulus squared
    // double real_part = creal(z);
    // double imag_part = cimag(z);
    // double mod_squared = real_part * real_part + imag_part * imag_part;

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
    const uint8_t pol, const double r, const double n_i, const double th_i
)
{
    double power = 0.0;

    // Placeholder, must update with proper computation!!!
    //
    // if (pol == 0)  // s-polarization
    // {
    //     power = (
    //         ( n_i * cos(th_i) - n_f * cos(th_f) )
    //         / ( n_i * cos(th_i) + n_f * cos(th_f) )
    //     );
    // } else  // p-polarization
    // {
    //     power = (
    //         ( n_f * cos(th_i) - n_i * cos(th_f) )
    //         / ( n_f * cos(th_i) + n_i * cos(th_f) )
    //     );
    // }

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
    const uint8_t polarization, const double n_i, const double n_f, const double th_i, const double th_f
)
{
    const double r = interface_r(polarization, n_i, n_f, th_i, th_f);

    // double R = R_from_r(r);
    // return R;

    return r;
}


/**
 * @brief
 *
 * @param
 * @param
 * @return
 */
double interface_T(
    const uint8_t polarization, const double n_i, const double n_f, const double th_i, const double th_f
)
{
    // const double t = interface_t(polarization, n_i, n_f, th_i, th_f);
    //
    // // double T = T_from_t(t);
    // // return T;
    //
    // return t;

    return 0.0;
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
    uint8_t pol, double n_list[], double d_list[], uint8_t num_layers, double th_0, double lam_vac
)
{
    // is_forward_angle()
    // list_snell()
    // interface_t()
    // interface_r()
    // make_2x2_array()
    // R_from_r()
    // T_from_t()
    // power_entering_from_r()
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
    uint8_t pol, double n_list[], double d_list[], double th_0, double lam_vac
)
{
    return 0.0;
}


/**
 * @brief Calculates ellipsometric parameters, in radians.
 *
 * Warning: Conventions differ. You may need to subtract pi/2 or whatever.
 *
 * @param n_list
 * @param d_list
 * @param th_0
 * @param lam_vac
 * @return
 */
double ellips(
    double n_list[], double d_list[], double th_0, double lam_vac
)
{
    return 0.0;
}


/**
 * @brief Calculates reflected and transmitted power for unpolarized light.
 *
 * @param n_list
 * @param d_list
 * @param th_0
 * @param lam_vac
 * @return
 */
double unpolarized_RT(
    double n_list[], double d_list[], double th_0, double lam_vac
)
{
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
 * @brief
 *
 * @param d_list The list of thicknesses of layers
 * @param d_list_size
 * @param final_answer
 * @return
 */
uint8_t layer_starts(
    const double d_list[], uint8_t d_list_size, double final_answer[]
)
{
    return 0;
}


/**
 * @brief An array listing what proportion of light is absorbed in each layer.
 *
 *
 *
 * @param coh_tmm_data
 * @return
 */
double absorp_in_each_layer(double coh_tmm_data)
{
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
