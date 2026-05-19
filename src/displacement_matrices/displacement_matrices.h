#ifndef DISPLACEMENT_MATRICES_H
#define DISPLACEMENT_MATRICES_H

#include "flint/gr.h"

// Displacement type, type
typedef enum { DISP_PLUS, DISP_MINUS } disp_type_t;

/**
 *
 * @brief This is the safe-book version relying purely on the optimisation of FLINT calculations
 * @deprecated This naive mathematical implementation is deprecated due to its
 * higher complexity. Please use `gr_mat_displacement()` instead, which
 * utilizes a more optimized O(nm) pointer-arithmetic approach. (Currently O(n^3) naive)
 * @param[out] D Resulting L shaped matrix of xnx
 * @param[in] A Input square nxn toeplitz matrix
 * @param[in] ctx The FLINT context object
 */
int gr_mat_displacement_square_safe(gr_mat_t D, gr_mat_t A, gr_ctx_t ctx);

/**
 * @brief Returns the displacement matrix nxm matrices.
 * Complexity: O(n*m), iterates through each element.
 * @param[out] D Resulting L shaped matrix
 * @param[in] A Input arbitrary matrix
 * @param[in] type Displacement type
 * @param[in] ctx The FLINT context object
 */
int gr_mat_displacement(gr_mat_t D, gr_mat_t A, disp_type_t type, gr_ctx_t ctx);

/**
 * @brief Returns the G and H matrices by the LU decomposition.
 * @warning G and H should be uninitialized.
 * @param[out] G Generator matrix G (left part of the displacement factorization)
 * @param[out] H Generator matrix H (right part of the displacement factorization)
 * @param[in] A Input Quasi Toeplitz or Toeplitz matrix
 * @param[in] type Displacement type
 * @param[in] ctx The FLINT context object
 */
int gr_mat_G_H(gr_mat_t G, gr_mat_t H, gr_mat_t A, disp_type_t type, gr_ctx_t ctx);

/**
 * @brief Reconstructs the original matrix A from its generators G and H.
 * Depending on the displacement type provided, it will automatically route
 * to either the Sigma LU (for DISP_PLUS) or Sigma UL (for DISP_MINUS)
 * mathematical reconstruction formula.
 * @param[out] A Pre-allocated dense matrix where the result will be stored
 * @param[in] G Generator matrix G (left part of the displacement factorization)
 * @param[in] H Generator matrix H (right part of the displacement factorization)
 * @param[in] type Displacement type
 * @param[in] ctx The FLINT context object
 */
int gr_mat_reconstruct_A(gr_mat_t A, gr_mat_t G, gr_mat_t H, disp_type_t type, gr_ctx_t ctx);
int gr_mat_reconstruct_A_v2(gr_mat_t A, gr_mat_t G, gr_mat_t H, disp_type_t type, gr_ctx_t ctx);

#endif
