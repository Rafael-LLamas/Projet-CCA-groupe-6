#include "addition.h"
#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "multiplication.h"

/*
Algorithm 10.1 Algorithme de type Strassen pour inverser une matrice
quasi-Toeplitz from [1].

Only square and non singular matrices are invertible. We will reject
any input of a non square matrix.
*/

int gr_toeplitz_inverse(gr_mat_t G_D, gr_mat_t H_D, gr_mat_t G_A, gr_mat_t H_A, gr_ctx_t ctx) {
  int res = GR_SUCCESS;
  // TODO: add a check of its determinant once implemented the toeplitz version (det must not be zero)

  // base case
  if (gr_mat_nrows(G_A, ctx) == 1) {
    // Gd[0, 0] = 1, Gd[0, 1] = Ga[0, 0] * Ha[0, 0]
    // Hd[0, 0] = 1, Hd[0, 1] = 0
    gr_one(gr_mat_entry(G_D, 0, 0, ctx), ctx);
    gr_mul(gr_mat_entry(G_D, 0, 1, ctx), gr_mat_entry(G_A, 0, 0, ctx), gr_mat_entry(H_A, 0, 0, ctx), ctx);
    gr_one(gr_mat_entry(H_D, 0, 0, ctx), ctx);
    gr_zero(gr_mat_entry(H_D, 0, 1, ctx), ctx);

    // if rank > 2 better to put everything else to 0
    for (slong i = 2; i < gr_mat_ncols(G_A, ctx); i++) {
      gr_zero(gr_mat_entry(G_D, 0, i, ctx), ctx);
      gr_zero(gr_mat_entry(H_D, 0, i, ctx), ctx);
    }
    return res;
  }

  return res;
}
