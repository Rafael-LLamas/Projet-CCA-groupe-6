cmake_minimum_required(VERSION 3.10.0)
project(model_project VERSION 0.1.0 LANGUAGES C)

include(CTest)
enable_testing()

add_compile_options(-g -O3 -Wall)
set(CMAKE_C_STANDARD 99)

# Setting path location
set(FLINT_SEARCH_PATHS
    "/usr/local"
    "/opt/homebrew"
    "/usr"
)

find_library(GMP_LIB gmp HINTS ${FLINT_SEARCH_PATHS} PATH_SUFFIXES lib)
find_path(GMP_INCLUDE gmp.h HINTS ${FLINT_SEARCH_PATHS} PATH_SUFFIXES include)

find_library(MPFR_LIB mpfr HINTS ${FLINT_SEARCH_PATHS} PATH_SUFFIXES lib)
find_path(MPFR_INCLUDE mpfr.h HINTS ${FLINT_SEARCH_PATHS} PATH_SUFFIXES include)

find_library(FLINT_LIB flint HINTS ${FLINT_SEARCH_PATHS} PATH_SUFFIXES lib)
find_path(FLINT_INCLUDE flint/flint.h HINTS ${FLINT_SEARCH_PATHS} PATH_SUFFIXES include)

#check
if(NOT GMP_LIB OR NOT FLINT_LIB OR NOT MPFR_LIB)
    message(FATAL_ERROR "Could not find one of the required libraries (GMP, MPFR, FLINT). Check your paths!!!")
endif()

# interface
add_library(flint_core INTERFACE)
target_include_directories(flint_core INTERFACE
    ${FLINT_INCLUDE}
    ${GMP_INCLUDE}
    ${MPFR_INCLUDE}
)
target_link_libraries(flint_core INTERFACE
    ${FLINT_LIB}
    ${GMP_LIB}
    ${MPFR_LIB}
)

####
# EXECUTABLES
###
set(PROJECT_SOURCES
    src/displacement_matrices/displacement_matrices.c
    src/random_toeplitz/random_toeplitz.c
    src/toeplitz_inverse/toeplitz_inverse.c
    src/utility/matrix_aux.c
    src/operations/multiplication.c
    src/operations/addition.c
)
add_library(project_lib STATIC ${PROJECT_SOURCES})
target_link_libraries(project_lib PUBLIC flint_core)
target_include_directories(project_lib PUBLIC
    src/displacement_matrices
    src/random_toeplitz
    src/toeplitz_inverse
    src/utility
    src/operations
)
add_executable(main src/main.c)
target_link_libraries(main PRIVATE project_lib)

####
# Tests
###
add_executable(test_displacement tests/displacement.test.c)
target_link_libraries(test_displacement PRIVATE project_lib)
add_test(NAME Displacement.toeplitz_deplacement                     COMMAND test_displacement toeplitz_deplacement)
add_test(NAME Displacement.manual_deplacement                       COMMAND test_displacement manual_deplacement)
add_test(NAME Displacement.large_deplacement                        COMMAND test_displacement large_deplacement)
add_test(NAME Displacement.toeplitz_G_H                             COMMAND test_displacement toeplitz_G_H)
add_test(NAME Displacement.quasi_toeplitz_G_H                       COMMAND test_displacement quasi_toeplitz_G_H)
add_test(NAME Displacement.toeplitz_reconstruction                  COMMAND test_displacement toeplitz_reconstruction)
add_test(NAME Displacement.quasi_toeplitz_reconstruction            COMMAND test_displacement quasi_toeplitz_reconstruction)

add_executable(test_utility tests/utility.test.c)
target_link_libraries(test_utility PRIVATE project_lib)
add_test(NAME Utility.LU_detatch                                    COMMAND test_utility LU_detatch)
add_test(NAME Utility.LU_detatch2                                   COMMAND test_utility LU_detatch2)

add_executable(test_random_toeplitz tests/random_toeplitz.test.c)
target_link_libraries(test_random_toeplitz PRIVATE project_lib)
add_test(NAME Random.gr_mat_random_toeplitz                                COMMAND test_random_toeplitz gr_mat_random_toeplitz)

add_executable(test_operations_multiplication tests/multiplication.test.c)
target_link_libraries(test_operations_multiplication PRIVATE project_lib)
add_test(NAME Operations.multiplication_generators                  COMMAND test_operations_multiplication multiplication_generators)

add_executable(test_operations_addition tests/addition.test.c)
target_link_libraries(test_operations_addition PRIVATE project_lib)
add_test(NAME Operations.addition_generators_toepltiz               COMMAND test_operations_addition addition_generators_toeplitz)

add_executable(test_inversion tests/inverse.test.c)
target_link_libraries(test_inversion PRIVATE project_lib)
add_test(NAME Inverse.base_case                                     COMMAND test_inversion inverse_base_case)
add_test(NAME Inverse.2x2                                           COMMAND test_inversion inverse_2x2)
