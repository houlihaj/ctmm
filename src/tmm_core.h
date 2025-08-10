//
// Created by johnh on 8/2/2025.
//

#ifndef TMM_CORE_H
#define TMM_CORE_H


#include <stdbool.h>


typedef struct {
    double complex poyn;
    double complex absor;
    double complex Ex;
    double complex Ey;
    double complex Ez;
} PositionResolvedData;


typedef struct {
    double psi;  // units of radians
    double delta;  // units of radians
} EllipsData;


uint8_t is_forward_angle(double complex n, double theta, bool* answer);


uint8_t snell(double n_1, double complex n_2, double th_1, double* th_2_guess);


uint8_t list_snell(
    const double complex n_list[],
    uint8_t n_list_size,
    double th_0,
    double complex angles[]
);


uint8_t interface_r(
    uint8_t polarization,
    double complex n_i,
    double complex n_f,
    double complex th_i,
    double complex th_f,
    double complex* r
);


uint8_t interface_t(
    uint8_t polarization,
    double complex n_i,
    double complex n_f,
    double complex th_i,
    double complex th_f,
    double complex* t
);


uint8_t R_from_r(double complex r, double* R);


uint8_t T_from_t(
    uint8_t pol,
    double complex t,
    double complex n_i,
    double complex n_f,
    double complex th_i,
    double complex th_f,
    double* T
);


uint8_t power_entering_from_r(
    uint8_t pol,
    double complex r,
    double complex n_i,
    double th_i,
    double* power
);


uint8_t interface_R(
    uint8_t polarization,
    double complex n_i,
    double complex n_f,
    double th_i,
    double th_f,
    double* R
);


uint8_t interface_T(
    uint8_t polarization,
    double complex n_i,
    double complex n_f,
    double th_i,
    double th_f,
    double* T
);


uint8_t coh_tmm(
    uint8_t pol,
    const double complex n_list[],
    const double d_list[],
    uint8_t num_layers,
    double th_0,
    double lam_vac,
    CohTmmData* coh_tmm_data
);


uint8_t coh_tmm_reverse(
    uint8_t pol,
    const double complex n_list[],
    const double d_list[],
    uint8_t num_layers,
    double th_0,
    double lam_vac,
    CohTmmData* coh_tmm_data
);


uint8_t ellips(
    const double complex n_list[],
    const double d_list[],
    uint8_t num_layers,
    double th_0,
    double lam_vac,
    EllipsData* ellips_data
);


uint8_t unpolarized_RT(
    const double complex n_list[],
    const double d_list[],
    uint8_t num_layers,
    double th_0,
    double lam_vac,
    double* R,
    double* T
);


uint8_t position_resolved(
    uint8_t layer,
    double distance,
    CohTmmData* coh_tmm_data,
    PositionResolvedData* position_resolved_data
);


uint8_t find_in_structure(
    const double d_list[],
    uint8_t d_list_size,
    double distance,
    double interface_info[]
);


uint8_t find_in_structure_with_inf(
    const double d_list[],
    uint8_t d_list_size,
    double distance,
    double interface_info[]
);


uint8_t layer_starts(
    const double d_list[], uint8_t d_list_size, double final_answer[]
);


uint8_t absorp_in_each_layer(
    CohTmmData* coh_tmm_data, uint8_t num_layers, double final_answer[]
);


// uint8_t inc_group_layers(double n_list[], double d_list[], double c_list[]);
//
//
// uint8_t inc_tmm(
//     uint8_t pol,
//     double complex n_list[],
//     double d_list[],
//     uint8_t c_list[],
//     uint8_t num_layers,
//     double th_0,
//     double lam_vac
// );
//
//
// uint8_t inc_absorp_in_each_layer(double inc_data);
//
//
// uint8_t inc_find_absorp_analytic_fn(uint8_t layer, double inc_data);


#endif //TMM_CORE_H
