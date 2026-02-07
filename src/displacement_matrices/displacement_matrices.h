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
 * @param rank Rank of the displacement matrix of A. (Attention for G and H the function will change: "-1" if already init., "0" if not init. )
 */
int gr_mat_G_H(gr_mat_t G, gr_mat_t H, gr_mat_t A, slong *rank, gr_ctx_t ctx);

#endif