cmake_minimum_required(VERSION 3.10.0)
project(model_project VERSION 0.1.0 LANGUAGES C)

include(CTest)
enable_testing()

add_compile_options(-g -O0 -Wall)
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

#find all files
file(GLOB_RECURSE PROJECT_SOURCES CONFIGURE_DEPENDS
    "src/gr_struct_mat/*.c"
)

#filter
list(FILTER PROJECT_SOURCES EXCLUDE REGEX "src/gr_struct_mat/test/.*")
list(FILTER PROJECT_SOURCES EXCLUDE REGEX "src/archive/.*")

add_library(project_lib STATIC ${PROJECT_SOURCES})

target_include_directories(project_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gr_struct_mat/utility
)
target_compile_options(project_lib PUBLIC
    "-I${FLINT_INCLUDE}"
    "-I${GMP_INCLUDE}"
    "-I${MPFR_INCLUDE}"
)

target_link_libraries(project_lib PUBLIC
    ${FLINT_LIB}
    ${GMP_LIB}
    ${MPFR_LIB}
)
add_executable(example src/example.c)
target_link_libraries(example PRIVATE project_lib)

####
# Tests
###

function(add_project_test TEST_TARGET FILE_NAME)
  add_executable(${TEST_TARGET} "src/gr_struct_mat/test/${FILE_NAME}")
  target_link_libraries(${TEST_TARGET} PRIVATE project_lib)
endfunction()

# reconstruction
add_project_test(test_reconstruction reconstruction.test.c)
add_test(NAME Reconstruction.toeplitz_reconstruction       COMMAND test_reconstruction toeplitz_reconstruction)
add_test(NAME Reconstruction.quasi_toeplitz_reconstruction COMMAND test_reconstruction quasi_toeplitz_reconstruction)

# Compression
add_project_test(test_compression compress.test.c)
add_test(NAME Compression.compress_addition_toeplitz       COMMAND test_compression compress_addition_toeplitz)
add_test(NAME Compression.compress_multiplication_toeplitz COMMAND test_compression compress_multiplication_toeplitz)
