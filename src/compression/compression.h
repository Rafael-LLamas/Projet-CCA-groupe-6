#ifndef COMPRESSION_H
#define COMPRESSION_H

#include "flint/gr.h"

/**
 * @brief the function will compress generators of A to D
 * @param[out] G_d Input the generator G of the matrix D
 * @param[out] H_d Input the generator H of the matrix D
 * @param[in] G_a Input the generator G of the matrix A
 * @param[in] H_a Input the generator H of the matrix A
 * @param[in] ctx The FLINT context object
 */
int gr_mat_generator_compress(gr_mat_t G_d, gr_mat_t H_d, gr_mat_t G_a, gr_mat_t H_a, gr_ctx_t ctx);

#endif
