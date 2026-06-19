#include <stdio.h>
#include <stdlib.h>

#include "../../gr_struct_mat.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"
#include "utility.h"

// Takes the dense matrix A and returns the generators G and H.
int _gr_mat_G_H(gr_mat_t G, gr_mat_t H, gr_mat_t A, disp_type_t type, gr_ctx_t ctx)
{

  /*
   * Returning the generators of a matrix A,
   * also can be called the rank optimization part of the toeplitz matrices.
   *
   * The displacement matrix of A (<- D) can be factored into G and H^T.
   * The product of G and H^T will give us the displacement matrix back.
   *
   * The mathematics definitions two types of displacement operators:
   *            Phi+ (DISP_PLUS):   D = A - Z * A * Z^T
   *            Phi- (DISP_MINUS):  D = A - Z^T * A * Z
   *
   * Regardless of which operator generated D, the resulting displacement matrix
   * can be expressed as a sum of outer products: D = G * H^T.
   *
   * Because our goal here is strictly to factorize D into a left part (G)
   * and a right part (H^T), we use an LU decomposition
   *
   * (D = P * L * U) for both types without changing the logic. By setting
   * G = P * L and H^T = U, we get G * H^T = D.
   */

  int status = GR_SUCCESS;
  slong m = gr_mat_nrows(A, ctx);
  slong n = gr_mat_ncols(A, ctx);

  slong *P = flint_malloc(m * sizeof(slong)); // Permutation table

  gr_mat_t D;
  gr_mat_init(D, m, n, ctx);

  status |= _gr_mat_displacement(D, A, type, ctx); // D <- displacement matrix

  if (status != GR_SUCCESS)
  {
    gr_mat_clear(D, ctx);
    flint_free(P);
    return status;
  }

  gr_mat_t LU, L, U;
  slong rank;

  gr_mat_init(LU, m, n, ctx);
  gr_mat_init(L, m, m, ctx);
  gr_mat_init(U, m, n, ctx);

  status |= gr_mat_lu(&rank, P, LU, D, 0, ctx); // LU decomposition

  if (rank < 1)
  {
    gr_mat_init(G, m, 1, ctx);
    gr_mat_init(H, n, 1, ctx);
    status |= gr_mat_zero(G, ctx);
    status |= gr_mat_zero(H, ctx);
    flint_free(P);
    gr_mat_clear(D, ctx);
    gr_mat_clear(LU, ctx);
    gr_mat_clear(L, ctx);
    gr_mat_clear(U, ctx);
    return status;
  }

  gr_mat_init(G, m, rank, ctx);
  gr_mat_init(H, n, rank, ctx);

  if (status == GR_SUCCESS)
  {
    status |= _gr_mat_lu_detach(L, U, LU, ctx); // detach the LU format to L and U

    // copy the values from L to G with correct permutation
    for (slong i = 0; i < m; i++)
      for (slong j = 0; j < rank; j++) status |= gr_set(gr_mat_entry_ptr(G, P[i], j, ctx), gr_mat_entry_srcptr(L, i, j, ctx), ctx);

    // copy the values from U to H
    for (slong i = 0; i < rank; i++)
      for (slong j = 0; j < n; j++) status |= gr_set(gr_mat_entry_ptr(H, j, i, ctx), gr_mat_entry_srcptr(U, i, j, ctx), ctx);
  }

  flint_free(P);
  gr_mat_clear(D, ctx);
  gr_mat_clear(LU, ctx);
  gr_mat_clear(L, ctx);
  gr_mat_clear(U, ctx);
  return status;
}