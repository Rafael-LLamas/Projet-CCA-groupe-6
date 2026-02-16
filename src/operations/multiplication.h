#ifndef MULTIPLICATION_H
#define MULTIPLICATION_H

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"

int test_multiplication_generateurs();

int gr_multiplication_generateur_deplacement(gr_mat_t G_c, gr_mat_t H_c, gr_mat_t G_a, gr_mat_t H_a, gr_mat_t G_b,
                                             gr_mat_t H_b, gr_ctx_t ctx);

int gr_multiplication_toeplitz(gr_mat_t C, gr_mat_t A, gr_mat_t B, gr_ctx_t ctx);

#endif