#ifndef TOEPLITZ_INVERSE_H
#define TOEPLITZ_INVERSE_H

#include "flint/gr.h"

/**
 * @brief Strassen algorithm to inverse a Toeplitz matrix.
 * @warning Matrix A must be square and non singular and G_D and H_D must be initialized *for now*.
 * @param[out] D Destination matrix represented by its G and H.
 * @param[in] A Input Toeplitz matrix represented by its G and H.
 */
int gr_toeplitz_inverse(gr_mat_t G_D, gr_mat_t H_D, gr_mat_t G_A, gr_mat_t H_A, gr_ctx_t ctx);

#endif
