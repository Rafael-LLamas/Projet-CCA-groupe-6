#include <stdio.h>
#include <stdlib.h>

#include "../gr_struct_mat.h"
#include "utility/utility.h"

#include "flint/flint.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"

int _gr_struct_mat_generator_compress(gr_mat_t G_D, gr_mat_t H_D, gr_mat_t G_A, gr_mat_t H_A, gr_ctx_t ctx)
{
  int status = GR_SUCCESS;

  slong n_G = gr_mat_nrows(G_A, ctx);
  slong n_H = gr_mat_nrows(H_A, ctx);
  slong k = gr_mat_ncols(G_A, ctx);

  gr_mat_t G_AT;
  gr_mat_init(G_AT, k, n_G, ctx);
  status |= gr_mat_transpose(G_AT, G_A, ctx);

  slong *P = flint_malloc(k * sizeof(slong)); // Permutation table
  slong rank;

  gr_mat_t LU, L, U;
  gr_mat_init(LU, k, n_G, ctx);
  gr_mat_init(L, k, k, ctx);
  gr_mat_init(U, k, n_G, ctx);

  // LU decomposition: P * GT = L * U
  status |= gr_mat_lu(&rank, P, LU, G_AT, 0, ctx);

  if (rank < 1 || status != GR_SUCCESS)
  {
    gr_mat_init(G_D, n_G, 1, ctx);
    gr_mat_init(H_D, n_H, 1, ctx);
    status |= gr_mat_zero(G_D, ctx);
    status |= gr_mat_zero(H_D, ctx);

    flint_free(P);
    gr_mat_clear(G_AT, ctx);
    gr_mat_clear(LU, ctx);
    gr_mat_clear(L, ctx);
    gr_mat_clear(U, ctx);
    return status;
  }

  status |= gr_mat_lu_detach(L, U, LU, ctx);

  gr_mat_init(G_D, n_G, rank, ctx);

  gr_mat_t H_perm, L_part;
  gr_mat_init(H_perm, n_H, k, ctx);
  gr_mat_init(L_part, k, rank, ctx);

  if (status == GR_SUCCESS)
  {
    // G_D = U^T for [0:r-1],*
    for (slong i = 0; i < rank; i++)
    {
      for (slong j = 0; j < n_G; j++) { status |= gr_set(gr_mat_entry_ptr(G_D, j, i, ctx), gr_mat_entry_srcptr(U, i, j, ctx), ctx); }
    }

    // permutation to H into H_perm
    for (slong i = 0; i < k; i++)
    {
      for (slong r = 0; r < n_H; r++)
      {
        status |= gr_set(gr_mat_entry_ptr(H_perm, r, i, ctx), gr_mat_entry_srcptr(H_A, r, P[i], ctx), ctx);
      }
    }

    // truncate L to L_part ∗,[0:r−1]
    for (slong i = 0; i < k; i++)
    {
      for (slong j = 0; j < rank; j++) { status |= gr_set(gr_mat_entry_ptr(L_part, i, j, ctx), gr_mat_entry_srcptr(L, i, j, ctx), ctx); }
    }

    // calculate H_D = H_perm * L_part
    gr_mat_init(H_D, n_H, rank, ctx);
    status |= gr_mat_mul(H_D, H_perm, L_part, ctx);
  }

  flint_free(P);
  gr_mat_clear(G_AT, ctx);
  gr_mat_clear(LU, ctx);
  gr_mat_clear(L, ctx);
  gr_mat_clear(U, ctx);
  gr_mat_clear(H_perm, ctx);
  gr_mat_clear(L_part, ctx);

  return status;
}

int gr_struct_mat_compress(gr_struct_mat_t mat, gr_ctx_t ctx)
{
  int status = GR_SUCCESS;

  // Might upgrade the whole helper to use the structure instead

  gr_mat_t G1, H1;
  status |= _gr_struct_mat_generator_compress(G1, H1, mat->G, mat->H, ctx);

  gr_mat_t G2, H2;
  status |= _gr_struct_mat_generator_compress(H2, G2, mat->H, mat->G, ctx);

  gr_mat_swap(mat->G, G2, ctx);
  gr_mat_swap(mat->H, H2, ctx);

  gr_mat_clear(G1, ctx);
  gr_mat_clear(H1, ctx);
  gr_mat_clear(G2, ctx);
  gr_mat_clear(H2, ctx);

  return status;
}