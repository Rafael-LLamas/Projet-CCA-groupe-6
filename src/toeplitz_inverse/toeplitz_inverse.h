#ifndef TOEPLITZ_INVERSE_H
#define TOEPLITZ_INVERSE_H

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"

/**
 * @brief Strassen algorithm to inverse a Toeplitz matrix.
 * @param[out] D Destination matrix represented by its G and H.
 * @param[in] A Input Toeplitz matrix represented by its G and H.
 */
int gr_toeplitz_inverse(gr_mat_t G_D, gr_mat_t H_D, gr_mat_t G_A, gr_mat_t H_D, gr_ctx_t ctx);

#endif

