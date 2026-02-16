#ifndef ADDITION_H
#define ADDITION_H

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"

int gr_mat_addition_generateur(gr_mat_t G_a, gr_mat_t H_a, gr_mat_t G_b, gr_mat_t H_b, gr_mat_t G_c, gr_mat_t H_c,
                               gr_ctx_t ctx);

#endif