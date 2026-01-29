#ifndef RANDOM_TOEPLITZ_H
#define RANDOM_TOEPLITZ_H

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"

int test_random_toeplitz();

void random_toeplitz(gr_mat_t A, int n, int m, flint_rand_t state, gr_ctx_t ctx);

#endif