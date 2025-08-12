#include <stdio.h>
#include <stdlib.h>


int main(void) {
    const char *filename = "data.csv";  // replace with your actual file name
    FILE *file = fopen(filename, "r");

    if (file == NULL) {
        perror("Error opening file");
        return EXIT_FAILURE;
    }

    char buffer[1024]; // Buffer to hold each line (adjust size if needed)
    int line_count = 0;

    while (fgets(buffer, sizeof(buffer), file) != NULL) {
        line_count++;
    }

    fclose(file);

    printf("Total lines: %d\n", line_count);
    return EXIT_SUCCESS;
}
