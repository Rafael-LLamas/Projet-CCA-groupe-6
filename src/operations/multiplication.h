#ifndef MULTIPLICATION_H
#define MULTIPLICATION_H

#include "flint/gr.h"

int test_multiplication_generateurs();

int gr_multiplication_generateur_deplacement(gr_mat_t G_c, gr_mat_t H_c, gr_mat_t G_a, gr_mat_t H_a, gr_mat_t G_b,
                                             gr_mat_t H_b, gr_ctx_t ctx);

int gr_multiplication_generateur_deplacement_fast(gr_mat_t G_c, gr_mat_t H_c, gr_mat_t G_a, gr_mat_t H_a, gr_mat_t G_b,
                                                  gr_mat_t H_b, gr_ctx_t ctx);
int gr_mat_apply_struct_fast(gr_mat_t Res, gr_mat_t G, gr_mat_t H, gr_mat_t X, gr_ctx_t ctx);
#endif
