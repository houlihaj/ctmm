//
// Created by johnh on 8/2/2025.
//

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

#include "src/tmm_coherent.h"
#include "src/tmm_core.h"


#define MAX_LINE_LEN 1024
#define MAX_TOKENS 3  // 10


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

    printf("Total lines: %d\n", line_count);
    return 0;
}


int main(int argc, char** argv) {
    printf("Hello, World!\nMy name is John\n");

    ////////////////////////////////////////////////////////////////////////////
    // Read the complex refractive indices from disk
    FILE *fp = fopen(
        "../data/nk/gold-Johnson.csv",
        "r"
    );
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
    while (fgets(line, sizeof(line), fp)) {
        // Remove newline if present
        line[strcspn(line, "\r\n")] = 0;

        char *tokens[MAX_TOKENS];
        int i = 0;

        char *token = strtok(line, ",");
        while (token && i < MAX_TOKENS) {
            tokens[i++] = token;
            token = strtok(NULL, ",");
        }

        // Now you can use tokens[0], tokens[1], etc.
        if (i == 3) {
            float wavelength = atof(tokens[0]);
            float refractive_index = atof(tokens[1]);
            float extinction_coefficient = atof(tokens[2]);

            printf(
                "Wavelength: %.4f, Refractive Index (n): %.2f, Extinction Coefficient (k): %.3f\n",
                wavelength, refractive_index, extinction_coefficient
            );
        }

        count++;
    }

    printf("%d", count);

    fclose(fp);

    ////////////////////////////////////////////////////////////////////////////
    // Execute coh_tmm() to compute reflection and transmission coefficients

    // CohTmmData coh_tmm_data;  // allocate struct in outer scope
    // const uint8_t num_layers = 4;
    // CohTmmData_create(&coh_tmm_data, num_layers);
    //
    // // list of complex refractive indices
    // const double complex n_list[] = {
    //     1.0 + 0.0 * I, 2.2 + 0.0 * I, 3.3 + 0.3 * I, 1.0 + 0.0 * I
    // };
    //
    // // list of layer thicknesses in nm
    // const double d_list[] = {INFINITY, 100, 300, INFINITY};  // must start and end with INFINITY
    //
    // coh_tmm(0, n_list, d_list, num_layers, 0.0, 1550.0, &coh_tmm_data);
    //
    // printf("reflection coefficient, r = %.3f + %.3fi\n", creal(coh_tmm_data.r), cimag(coh_tmm_data.r));
    // printf("reflectivity, R = %.3f\n", coh_tmm_data.R);
    //
    // printf("transmission coefficient, t = %.3f + %.3fi\n", creal(coh_tmm_data.t), cimag(coh_tmm_data.t));
    // printf("transmissivity, T = %.3f\n", coh_tmm_data.T);
    //
    // CohTmmData_destroy(&coh_tmm_data);

    return 0;
}
