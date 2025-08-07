//
// Created by johnh on 8/2/2025.
//

#ifndef TMM_CORE_H
#define TMM_CORE_H


// uint8_t make_2x2_array(
//     const double a, const double b, const double c, const double d
// );
//
//
// uint8_t is_forward_angle(double complex n, double theta);
//
//
uint8_t snell(const double n_1, const double complex n_2, const double th_1, double* th_2_guess);
//
//
// uint8_t list_snell(
//     const double n_list[], const uint8_t n_list_size, const double th_0, double angles[]
// );
//
//
uint8_t interface_r(
    uint8_t polarization, double complex n_i, double complex n_f, double th_i, double th_f, double complex* r
);


uint8_t interface_t(
    uint8_t polarization, double complex n_i, double complex n_f, double th_i, double th_f, double complex* t
);


uint8_t R_from_r(const double complex r, double* R);
//
//
// double T_from_t(
//     const uint8_t pol, const double t, const double n_i, const double n_f, const double th_i, const double th_f
// );
//
//
// double power_entering_from_r(
//     const uint8_t pol, const double r, const double n_i, const double th_i
// );
//
//
// double interface_R(
//     const uint8_t polarization, const double n_i, const double n_f, const double th_i, const double th_f
// );
//
//
// double interface_T(
//     const uint8_t polarization, const double n_i, const double n_f, const double th_i, const double th_f
// );
//
//
// uint8_t coh_tmm(
//     uint8_t pol, double n_list[], double d_list[], uint8_t num_layers, double th_0, double lam_vac
// );
//
//
// double coh_tmm_reverse(
//     uint8_t pol, double n_list[], double d_list[], double th_0, double lam_vac
// );
//
//
// double ellips(
//     double n_list[], double d_list[], double th_0, double lam_vac
// );
//
//
// double unpolarized_RT(
//     double n_list[], double d_list[], double th_0, double lam_vac
// );
//
//
// double position_resolved(
//     uint8_t layer, double distance, double coh_tmm_data
// );
//
//
// uint8_t find_in_structure(
//     const double d_list[], const uint8_t d_list_size, const double distance
// );
//
//
// uint8_t find_in_structure_with_inf(
//     const double d_list[], const double distance
// );
//
//
// uint8_t layer_starts(
//     const double d_list[], uint8_t d_list_size, double final_answer[]
// );
//
//
// double absorp_in_each_layer(double coh_tmm_data);
//
//
// double inc_group_layers(double n_list[], double d_list[], double c_list[]);
//
//
// double inc_tmm(
//     uint8_t pol, double n_list[], double d_list[], double c_list[], double th_0, double lam_vac
// );
//
//
// double inc_absorp_in_each_layer(double inc_data);
//
//
// double inc_find_absorp_analytic_fn(uint8_t layer, double inc_data);


#endif //TMM_CORE_H
