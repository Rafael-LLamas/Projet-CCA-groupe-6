#include <stdio.h>
#include <stdlib.h>

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"
#include "gr_struct_mat.h"
#include "utility.h"

// Takes the dense matrix A and returns the generators G and H.
int _gr_mat_G_H(gr_mat_t G, gr_mat_t H, gr_mat_t A, structure_type_t struct_t, disp_type_t type, gr_ctx_t ctx)
{
  int status = GR_SUCCESS;
  slong m = gr_mat_nrows(A, ctx);
  slong n = gr_mat_ncols(A, ctx);

  slong *P = flint_malloc(m * sizeof(slong)); // Permutation table

  gr_mat_t D;
  gr_mat_init(D, m, n, ctx);

  status |= _gr_mat_displacement(D, A, struct_t, type, ctx);

  if (status != GR_SUCCESS)
  {
    gr_mat_clear(D, ctx);
    flint_free(P);
    return status;
  }

  gr_mat_t LU;
  slong rank;

  gr_mat_init(LU, m, n, ctx);
  status |= gr_mat_lu(&rank, P, LU, D, 0, ctx); // LU decomposition

  if (rank < 1)
  {
    gr_mat_init(G, m, 0, ctx);
    gr_mat_init(H, n, 0, ctx);

    flint_free(P);
    gr_mat_clear(D, ctx);
    gr_mat_clear(LU, ctx);
    return status;
  }

  gr_mat_init(G, m, rank, ctx);
  gr_mat_init(H, n, rank, ctx);

  if (status == GR_SUCCESS)
  {
    for (slong i = 0; i < m; i++)
      for (slong j = 0; j < rank; j++)
        if (i > j)
          status |= gr_set(gr_mat_entry_ptr(G, P[i], j, ctx), gr_mat_entry_srcptr(LU, i, j, ctx), ctx);
        else if (i == j)
          status |= gr_one(gr_mat_entry_ptr(G, P[i], j, ctx), ctx);
        else
          status |= gr_zero(gr_mat_entry_ptr(G, P[i], j, ctx), ctx);

    // U^T: H[j,i] = U[i,j]
    for (slong j = 0; j < n; j++)
      for (slong i = 0; i < rank; i++)
        if (i <= j)
          status |= gr_set(gr_mat_entry_ptr(H, j, i, ctx), gr_mat_entry_srcptr(LU, i, j, ctx), ctx);
        else
          status |= gr_zero(gr_mat_entry_ptr(H, j, i, ctx), ctx);
  }

  flint_free(P);
  gr_mat_clear(D, ctx);
  gr_mat_clear(LU, ctx);

  return status;
}