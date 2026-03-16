#ifndef DISPLACEMENT_MATRICES_H
#define DISPLACEMENT_MATRICES_H

#include "flint/gr.h"

// Displacement type, type
typedef enum { DISP_PLUS, DISP_MINUS } disp_type_t;

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
 * @param[in] type Displacement type (DISP_PLUS or DISP_MINUS)
 */
int gr_mat_displacement(gr_mat_t D, gr_mat_t A, disp_type_t type, gr_ctx_t ctx);

/**
 * @brief Returns the G and H matrices by the LU decomposition.
 * @warning G and H should be uninitialized.
 * @param[out] G Resulting G matrix
 * @param[out] H Resulting H^T matrix
 * @param[in] A Input Quasi Toeplitz or Toeplitz matrix
 */
int gr_mat_G_H(gr_mat_t G, gr_mat_t H, gr_mat_t A, disp_type_t type, gr_ctx_t ctx);

/**
 * @brief Returns the Matrix A from generators G and H via SigmaLU.
 * @param[in] G Generator M
 */
int gr_mat_reconstruct_A(gr_mat_t A, gr_mat_t G, gr_mat_t H, disp_type_t type, gr_ctx_t ctx);

#endif
