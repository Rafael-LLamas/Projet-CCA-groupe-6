#include <stdio.h>
#include <stdlib.h>

#include "flint/gr_types.h"

int gr_mat_generator_compress(gr_mat_t G_d, gr_mat_t H_d, gr_mat_t G_a, gr_mat_t H_a, gr_ctx_t ctx) {

  /*
   * This is an auxilary function to perform compression on Toeplitz generators
   *
   * Main goal is to achieve generators of D such that
   *        reconstruct(D) = reconstruct(A)
   * meanwhile having the rank of generators D less or equal to A.
   */

  int status = GR_UNABLE;
  return status;
}
