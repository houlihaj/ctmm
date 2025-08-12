#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>


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

    *num_lines = line_count;
    return EXIT_SUCCESS;
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
    float wavelengths[num_lines - 1];  // one less because of the header row
    float refractive_indices[num_lines - 1];
    float extinction_coefficients[num_lines - 1];

    ////////////////////////////////////////////////////////////////////////////
    // Read the complex refractive indices from disk
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
            wavelengths[i] = wavelength;
            float refractive_index = atof(tokens[1]);
            refractive_indices[i] = refractive_index;
            float extinction_coefficient = atof(tokens[2]);
            extinction_coefficients[i] = extinction_coefficient;

            printf(
                "Wavelength: %.4f, Refractive Index (n): %.2f, Extinction Coefficient (k): %.3f\n",
                wavelength, refractive_index, extinction_coefficient
            );
        }

        count++;
    }

    printf("%d", count);

    fclose(fp);

    return 0;
}
