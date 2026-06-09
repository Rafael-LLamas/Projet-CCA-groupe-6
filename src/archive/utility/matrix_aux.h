#ifndef MATRIX_AUX_H
#define MATRIX_AUX_H

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"

#define GNMOD 1009

// to retire one day
#define FLINT_CHECK(x)                                                                                                 \
  do {                                                                                                                 \
    int _status = (x);                                                                                                 \
    if (_status != GR_SUCCESS) return _status;                                                                         \
  } while (0)

/**
 * @brief Fills an initialized matrix with random elements
 * Assigns the matrix D with random UINTS st < limit
 * @pre The matrix D must be previously initialized
 * @param[out] D The result matrix to be populated
 * @param[in] state The initialized FLINT random state
 * @param[in] limit The strictly upper bound for the random elements
 * @param[in] ctx The FLINT generic ring context
 */
int gr_mat_random_manual(gr_mat_t D, flint_rand_t state, ulong limit, gr_ctx_t ctx);

/**
 * @brief Detaches an LU decomposed matrix into its separate L and U components.
 * @pre L, U, and LU must be initialized. L and U must have dimensions compatible with LU.
 * @param[out] L The lower triangular result matrix.
 * @param[out] U The upper triangular result matrix.
 * @param[in] LU The combined LU matrix source.
 * @param[in] ctx The FLINT generic ring context.
 */
int gr_mat_lu_detach(gr_mat_t L, gr_mat_t U, gr_mat_t LU, gr_ctx_t ctx);

/**
 * @brief Applies the Z matrix (downward shift operator) on matrix M.
 * @pre Res and M must be initialized and have the same dimensions.
 * @param[out] Res Result matrix containing the shifted result.
 * @param[in] M Source to be shifted.
 * @param[in] ctx The FLINT generic ring context.
 */
int gr_mat_apply_Z(gr_mat_t Res, gr_mat_t M, gr_ctx_t ctx);

/**
 * @brief Applies the Z^T matrix (upwards shift operator) on matrix M.
 * @pre Res and M must be initialized and have the same dimensions.
 * @param[out] Res Result matrix containing the shifted result.
 * @param[in] M Source to be shifted.
 * @param[in] ctx The FLINT generic ring context.
 */
int gr_mat_apply_Zt(gr_mat_t Res, gr_mat_t M, gr_ctx_t ctx);

#endif