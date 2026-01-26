cmake_minimum_required(VERSION 3.10.0)
project(model_project VERSION 0.1.0 LANGUAGES C)

include(CTest)
enable_testing()

add_compile_options(-g -O0 -Wall)
set(CMAKE_C_STANDARD 99)

# Setting path location
set(FLINT_ROOT "/usr/local")
if(EXISTS "/opt/homebrew") # If apple silicon
    set(FLINT_ROOT "/opt/homebrew")
endif()

# Find GMP
find_library(GMP_LIB gmp HINTS ${FLINT_ROOT}/lib)
find_path(GMP_INCLUDE gmp.h HINTS ${FLINT_ROOT}/include)

# Find MPFR
find_library(MPFR_LIB mpfr HINTS ${FLINT_ROOT}/lib)
find_path(MPFR_INCLUDE mpfr.h HINTS ${FLINT_ROOT}/include)

# Find FLINT
find_library(FLINT_LIB flint HINTS ${FLINT_ROOT}/lib)
find_path(FLINT_INCLUDE flint/flint.h HINTS ${FLINT_ROOT}/include)

add_executable(test_bordel 
    "tests/Test_Bordel_Flint.c"
)
add_executable(random_toeplitz 
    "src/random_toeplitz/random_toeplitz.c"
)

target_include_directories(test_bordel PRIVATE 
    ${FLINT_INCLUDE} 
    ${GMP_INCLUDE} 
    ${MPFR_INCLUDE}
)

target_link_libraries(test_bordel PRIVATE 
    ${FLINT_LIB}
    ${GMP_LIB}
    ${MPFR_LIB}
)

target_include_directories(random_toeplitz PRIVATE 
    ${FLINT_INCLUDE} 
    ${GMP_INCLUDE} 
    ${MPFR_INCLUDE}
)

target_link_libraries(random_toeplitz PRIVATE 
    ${FLINT_LIB}
    ${GMP_LIB}
    ${MPFR_LIB}
)
