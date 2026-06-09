#include <stdlib.h>
#include <time.h>

#include "flint/flint.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"

int gr_mat_addition_generateur(gr_mat_t T, gr_mat_t U, gr_mat_t G, gr_mat_t H, gr_mat_t G_c, gr_mat_t H_c,
                               gr_ctx_t ctx) {

  int status = GR_SUCCESS;

  gr_mat_init(G_c, gr_mat_nrows(T, ctx), gr_mat_ncols(T, ctx) + gr_mat_ncols(G, ctx), ctx);
  gr_mat_init(H_c, gr_mat_nrows(U, ctx), gr_mat_ncols(U, ctx) + gr_mat_ncols(H, ctx), ctx);

  status |= gr_mat_concat_horizontal(G_c, T, G, ctx);
  status |= gr_mat_concat_horizontal(H_c, U, H, ctx);
  return status;
}
