#ifndef DISPLACEMENT_MATRICES_H
#define DISPLACEMENT_MATRICES_H

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"

#define FLINT_CHECK(x)                                                                                                 \
  do {                                                                                                                 \
    int _status = (x);                                                                                                 \
    if (_status != GR_SUCCESS) return _status;                                                                         \
  } while (0)

/**
 * @brief Main testing function for displacement matrices
 */
int test_displacement_matrices();

/**
 * @brief This is the safe-book version relying purely on the optimisation of FLINT calculations
 * Complexity: O(n^3) naive, but this is flint so probably less
 * @param[out] D Resulting L shaped matrix of xnx
 * @param[in] A Input square nxn toeplitz matrix
 */
int gr_mat_displacement_square(gr_mat_t D, gr_mat_t A, gr_ctx_t ctx);

/**
 * @brief Returns the displacement matrix nxm matrices.
 * Complexity: O(n*m), iterates through each element.
 * @param[out] D Resulting L shaped matrix
 * @param[in] A Input matrix
 */
int gr_mat_displacement(gr_mat_t D, gr_mat_t A, gr_ctx_t ctx);

#endif