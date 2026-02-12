cmake_minimum_required(VERSION 3.10.0)
project(model_project VERSION 0.1.0 LANGUAGES C)

include(CTest)
enable_testing()

add_compile_options(-g -O0 -Wall)
set(CMAKE_C_STANDARD 99)

# Setting path location
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
set(SOURCES
    src/main.c
    src/displacement_matrices/displacement_matrices.c
    src/random_toeplitz/random_toeplitz.c
    src/utility/matrix_aux.c
    src/multiplication_toeplitz/multiplication.c
)

add_executable(main ${SOURCES})
target_link_libraries(main PRIVATE flint_core)
target_include_directories(main PRIVATE # this one helps with include paths
    src/displacement_matrices
    src/random_toeplitz
    src/utility
    src/multiplication
)

####
# Tests
###

add_test(NAME DisplacementTest      COMMAND main displacement)
add_test(NAME RandomTest            COMMAND main random)
add_test(NAME AuxTest               COMMAND main aux)
add_test(NAME MultiplicationTest    COMMAND main multiplication)