#include <stdlib.h>
#include <time.h>

#include "flint/flint.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"

int gr_mat_addition_generateur(gr_mat_t G_a, gr_mat_t H_a, gr_mat_t G_b, gr_mat_t H_b, gr_mat_t G_c, gr_mat_t H_c,
                               gr_ctx_t ctx) {

  int status = GR_SUCCESS;

  gr_mat_init(G_c, gr_mat_nrows(G_a, ctx), gr_mat_ncols(G_a, ctx) + gr_mat_ncols(G_b, ctx), ctx);
  gr_mat_init(H_c, gr_mat_nrows(H_a, ctx), gr_mat_ncols(H_a, ctx) + gr_mat_ncols(H_b, ctx), ctx);

  status |= gr_mat_concat_horizontal(G_c, G_a, G_b, ctx);
  status |= gr_mat_concat_horizontal(H_c, H_a, H_b, ctx);
  return status;
}
