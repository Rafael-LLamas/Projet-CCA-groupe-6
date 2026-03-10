#ifndef RANDOM_TOEPLITZ_H
#define RANDOM_TOEPLITZ_H

#include "flint/gr_mat.h"

int random_quasi_toeplitz(gr_mat_t A, int n, int m, flint_rand_t state, gr_ctx_t ctx);
int random_toeplitz(gr_mat_t A, int n, int m, flint_rand_t state, gr_ctx_t ctx);
int rand_quasi_toeplitz(gr_mat_t A, int n, int m, int nb_rand, gr_ctx_t ctx);

#endif
