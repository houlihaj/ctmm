#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define MAX_LINE_LEN 1024
#define MAX_TOKENS 10


/**
 * @brief Example data.csv
 *
 * Name,Age,Score
 * Alice,23,89.5
 * Bob,30,92.0
 * Charlie,28,85.0
 *
 */
int main() {
    FILE *fp = fopen("data.csv", "r");
    if (!fp) {
        perror("Failed to open file");
        return 1;
    }

    char line[MAX_LINE_LEN];

    // Read header line (optional)
    if (fgets(line, sizeof(line), fp)) {
        printf("Header: %s", line);
    }

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
            char *name = tokens[0];
            int age = atoi(tokens[1]);
            float score = atof(tokens[2]);

            printf("Name: %s, Age: %d, Score: %.2f\n", name, age, score);
        }
    }

    fclose(fp);
    return 0;
}
