#ifndef COMPRESSION_H
#define COMPRESSION_H

#include "flint/gr.h"

// With arithmetic operations with generators the rank of the generators
// inflate artifitally.
// Compression of these generators allow to efficiently store the same
// matrix in algorithms.

/**
 * @brief the function will compress generators of A to D
 * @param[out] G_d Input/Output the generator G of the matrix D
 * @param[out] H_d Input/Output the generator H of the matrix D
 * @param[in] ctx The FLINT context object
 */
int gr_mat_generator_compress(gr_mat_t G_d, gr_mat_t H_d, gr_ctx_t ctx);

#endif
