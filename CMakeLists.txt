cmake_minimum_required(VERSION 3.31)
project(ctmm C)

set(CMAKE_C_STANDARD 99)

add_executable(${PROJECT_NAME}
        main.c
        src/tmm_core.c
        src/tmm_core.h
        src/tmm_absorp_fn.c
        src/tmm_absorp_fn.h
        src/tmm_coherent.c
        src/tmm_coherent.h
        src/tmm_math.c
        src/tmm_math.h
        src/tmm_util.c
        src/tmm_util.h
        src/tmm_incoherent.c
        src/tmm_incoherent.h
        src/tmm_logging.c
        src/tmm_logging.h
        external/log.c/log.c
)

# Link with the math library (-lm)
target_link_libraries(${PROJECT_NAME} m)
