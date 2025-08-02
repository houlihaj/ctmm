cmake_minimum_required(VERSION 3.31)
project(ctmm C)

set(CMAKE_C_STANDARD 99)

add_executable(${PROJECT_NAME}
        main.c
        tmm_core.c
        tmm_core.h
        tmm_absorp_fn.c
        tmm_absorp_fn.h)

# Link with the math library (-lm)
target_link_libraries(${PROJECT_NAME} m)
