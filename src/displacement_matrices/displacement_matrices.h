#ifndef DISPLACEMENT_MATRICES_H
#define DISPLACEMENT_MATRICES_H

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"

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
int gr_mat_displacement_square_safe(gr_mat_t D, gr_mat_t A, gr_ctx_t ctx);

/**
 * @brief Returns the displacement matrix nxm matrices.
 * Complexity: O(n*m), iterates through each element.
 * @param[out] D Resulting L shaped matrix
 * @param[in] A Input matrix
 */
int gr_mat_displacement(gr_mat_t D, gr_mat_t A, gr_ctx_t ctx);

/**
 * @brief Returns the G and H matrices by the LU decomposition.
 * @param[out] G Resulting G matrix
 * @param[out] H Resulting H^T matrix
 * @param[in] A Input Quasi Toeplitz or Toeplitz matrix
 * @param rank There are 3 cases: (1-) If know the rank of A, the rank of A. (2-) If not know the rank and G and H are
 * initialized, -1. (3-) If not know the rank of A and G and H are uninitialized, the function will initialize, put 0.
 */
int gr_mat_G_H(gr_mat_t G, gr_mat_t H, gr_mat_t A, gr_ctx_t ctx);

/**
 * @brief Returns the Matrix A from generators G and H via SigmaLU.
 * @param[in] G Generator M
 */
int gr_mat_reconstruct_A_safe(gr_mat_t A, gr_mat_t G, gr_mat_t H, gr_ctx_t ctx);

#endif