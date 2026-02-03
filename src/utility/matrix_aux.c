#include "matrix_aux.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"

#include <stdlib.h>
#include <time.h>

int gr_mat_random_manual(gr_mat_t D, gr_ctx_t ctx) { return gr_not_implemented(); }

int gr_mat_lu_detach(gr_mat_t L, gr_mat_t U, gr_mat_t LU, gr_ctx_t ctx) { return gr_not_implemented(); }

// TODO - might integrate with cmake tests later
int test_matrix_aux() {
  flint_printf("*----------* Matrix Auxilary Test *----------*\n");
  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, 1009);
  int res = GR_SUCCESS;
  return res;
}