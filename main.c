#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <float.h>
#include <math.h>
#include "src/tmm_core.h"


#define SPEED_OF_LIGHT 299792458

#define MAX_LINE_LEN 1024
#define MAX_TOKENS 3  // 10  // number of data columns


/**
 * @brief Get the number of lines from the datafile
 *
 * @param filename
 * @param num_lines
 * @return
 */
int get_number_of_lines_in_datafile(char* filename, uint16_t* num_lines)
{
    FILE *file = fopen(filename, "r");

    if (file == NULL) {
        perror("Error opening file");
        return EXIT_FAILURE;
    }

    char buffer[MAX_LINE_LEN];  // buffer to hold each line (adjust size if needed)
    int line_count = 0;

    while (fgets(buffer, sizeof(buffer), file) != NULL) {
        line_count++;
    }

    fclose(file);

    // printf("Total lines: %d\n", line_count);

    *num_lines = line_count;
    return EXIT_SUCCESS;
}


/**
 * @brief Read the complex refractive indices from disk
 *
 * @param filename
 * @param num_lines
 * @param lamba_arr Wavelengths
 * @param n_arr Refractive indices
 * @param k_arr Extinction coefficients
 * @return
 */
int read_nk_data(
    char* filename,
    const uint16_t num_lines,
    float lamba_arr[],
    float n_arr[],
    float k_arr[]
)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Failed to open file");
        return 1;
    }

    char line[MAX_LINE_LEN];

    // Read header line (optional)
    if (fgets(line, sizeof(line), fp)) {
        printf("Header: %s", line);
    }

    int count = 1;  // including the header row in our count

    // Read data lines
    for (int i = 1; i < num_lines; i++)
    {
        fgets(line, sizeof(line), fp);

        // Remove newline if present
        line[strcspn(line, "\r\n")] = 0;

        char *tokens[MAX_TOKENS];
        int j = 0;

        char *token = strtok(line, ",");
        while (token && j < MAX_TOKENS) {
            tokens[j++] = token;
            token = strtok(NULL, ",");
        }

        // Now you can use tokens[0], tokens[1], etc.
        if (j == MAX_TOKENS) {
            float wavelength = atof(tokens[0]);
            lamba_arr[i] = wavelength;
            float refractive_index = atof(tokens[1]);
            n_arr[i] = refractive_index;
            float extinction_coefficient = atof(tokens[2]);
            k_arr[i] = extinction_coefficient;

            // printf(
            //     "Wavelength: %.4f, Refractive Index (n): %.2f, Extinction Coefficient (k): %.3f\n",
            //     wavelength, refractive_index, extinction_coefficient
            // );
        }

        count++;
    }

    // printf("%d", count);

    fclose(fp);

    return 0;
}


int main(void) {
    printf("Hello, World!\nMy name is John\n");

    // Declare the datafile path
    char* filename = "../data/nk/gold-Johnson.csv";

    ////////////////////////////////////////////////////////////////////////////
    // Get the number of lines from the datafile
    uint16_t num_lines;
    get_number_of_lines_in_datafile(filename, &num_lines);
    printf("num_lines == %d\n", num_lines);

    ////////////////////////////////////////////////////////////////////////////
    // Declare the arrays that will hold the data
    const uint16_t num_rows = num_lines - 1;  // one less because of the header row
    float lambda_list[num_rows];  // units of micron
    float n_dat[num_rows];
    float k_dat[num_rows];

    ////////////////////////////////////////////////////////////////////////////
    // Read the complex refractive indices from disk
    read_nk_data(
        filename, num_lines, lambda_list, n_dat, k_dat
    );

    ////////////////////////////////////////////////////////////////////////////
    // Complex refractive index (eta) array
    double complex eta_dat[num_rows];

    for (int i = 0; i < num_rows; i++)
    {
        eta_dat[i] = n_dat[i] + k_dat[i] * I;
    }

    const double d_list[] = {INFINITY, 0.300, INFINITY};  // units of micron

    double R_list[num_rows];
    double T_list[num_rows];

    CohTmmData coh_tmm_data;  // allocate struct in outer scope
    const uint8_t num_layers = 3;
    CohTmmData_create(&coh_tmm_data, num_layers);
    
    for (int i = 0; i < num_rows; i++)
    {
        double complex n_list[] = {
            1.0 + 0.0 * I, eta_dat[i], 1.5 + 0.0 * I
        };

        coh_tmm(
            0,  // const uint8_t pol,
            n_list,  // const double complex n_list[],
            d_list,  // const double d_list[],
            num_layers,  // const uint8_t num_layers,
            0,  // const double th_0,
            lambda_list[i],  // const double lam_vac,
            &coh_tmm_data  // CohTmmData* coh_tmm_data
        );

        R_list[i] = coh_tmm_data.R;
        T_list[i] = coh_tmm_data.T;
    }

    CohTmmData_destroy(&coh_tmm_data);

    for (int i = 0; i < num_rows; i++)
    {
        printf("reflectivity, R = %.3f\n", R_list[i]);
        // printf("transmissivity, T = %.3f\n", T_list[i]);
    }

    return 0;
}
