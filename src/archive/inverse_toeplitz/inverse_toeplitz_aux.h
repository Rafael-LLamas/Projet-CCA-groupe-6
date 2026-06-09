#ifndef INVERSE_TOEPLITZ_AUX_H
#define INVERSE_TOEPLITZ_AUX_H

#include "flint/gr.h"

/**
 * @brief Slices the global generators into 4 quadrants and applies bleed correction.
 * This function extracts the generators and appends the necessary correction vectors to fix this bleed.
 */
int gr_mat_split_quadrants(gr_mat_t G_a, gr_mat_t H_a, gr_mat_t G_b, gr_mat_t H_b, gr_mat_t G_c, gr_mat_t H_c,
                           gr_mat_t G_d, gr_mat_t H_d, gr_mat_t G_A, gr_mat_t H_A, gr_ctx_t ctx);

/**
 * @brief Packs the local inverses (x, y, z, t) into global generators and fixes reverse bleed.
 * This function calculates the negative boundary vectors and adds them to the generators to 
 * cancel out the global shift.
 */
int gr_mat_pack_quadrants(gr_mat_t G_D, gr_mat_t H_D, gr_mat_t G_x, gr_mat_t H_x, gr_mat_t G_y, gr_mat_t H_y,
                          gr_mat_t G_z, gr_mat_t H_z, gr_mat_t G_t, gr_mat_t H_t, gr_ctx_t ctx);

#endif
