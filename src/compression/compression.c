#include <stdio.h>
#include <stdlib.h>

#include "displacement_matrices.h"

#include "flint/flint.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"
#include "flint/gr_types.h"

int gr_mat_generator_compress(gr_mat_t G_d, gr_mat_t H_d, gr_ctx_t ctx) {

  /*
   * This is an auxilary function to perform compression on Toeplitz generators
   *
   * Main goal is to achieve generators of D such that
   *        reconstruct(D) = reconstruct(A)
   * meanwhile having the rank of generators D less or equal to A.
   *
   * This code actually has no optimizations whatsoever and is only
   * a naive first implementation to get strassen up and working.
   */

  int status = GR_SUCCESS;
  slong n = gr_mat_nrows(G_d, ctx);
  slong m = gr_mat_nrows(H_d, ctx); // H is n x rank, so nrows = n too

  gr_mat_t A;
  gr_mat_init(A, n, m, ctx);

  status |= gr_mat_reconstruct_A(A, G_d, H_d, DISP_PLUS, ctx);
  if (status != GR_SUCCESS) {
    gr_mat_clear(A, ctx);
    return status;
  }

  gr_mat_clear(G_d, ctx);
  gr_mat_clear(H_d, ctx);

  status |= gr_mat_G_H(G_d, H_d, A, DISP_PLUS, ctx);

  gr_mat_clear(A, ctx);
  return status;
}