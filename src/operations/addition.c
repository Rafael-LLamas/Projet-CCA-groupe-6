#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "matrix_aux.h"
#include "random_toeplitz.h"

#include <stdlib.h>
#include <time.h>

int gr_mat_addition_generateur(gr_mat_t G_a, gr_mat_t H_a, gr_mat_t G_b, gr_mat_t H_b, gr_mat_t G_c, gr_mat_t H_c,
                               gr_ctx_t ctx) {

  slong new_rank = gr_mat_ncols(G_a, ctx) + gr_mat_ncols(G_b, ctx);
  gr_mat_init(H_c, gr_mat_nrows(G_a, ctx), new_rank, ctx);
  gr_mat_init(G_c, gr_mat_nrows(G_a, ctx), new_rank, ctx);
  /*
  C'est juste la concaténation des générateurs de déplacement de A et B
  Il n'y a rien de magique et complexe - dit rafael hamas
  */
  FLINT_CHECK(gr_mat_concat_horizontal(G_c, G_a, G_b, ctx));
  FLINT_CHECK(gr_mat_concat_horizontal(H_c, H_a, H_b, ctx));
  return GR_SUCCESS;
}