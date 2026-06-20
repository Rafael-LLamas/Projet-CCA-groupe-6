#include <stdio.h>
#include <stdlib.h>

#include "gr_struct_mat.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "flint/gr_types.h"


int _gr_mat_reconstruct(gr_mat_t A, gr_mat_t G, gr_mat_t H, disp_type_t type, gr_ctx_t ctx)
{

  /* Reconstructs the original A from its generators
   *
   * This function performs the inverse of the displacement operator
   * It is a computational approach to the multiplication of L and U using
   * diagonal pointer arithmetics
   *
   * The reconstruction depends on the displacement type used to generate
   * G and H:
   *
   * DISP_PLUS | Phi+ (Sigma LU)
   * Toeplitz data flows from the top-left corner, to compute cell A[i,j],
   * following the data backwards by walking diagonally UP and LEFT towards
   * the (0,0). We accumulate G[i-x] * H[j-x] and walk until we hit
   * the top or left wall (max_x = min(i, j)).
   *
   * DISP_MINUS | Phi- (Sigma UL)
   * A = Sum( U(rev(x_k)) * L(rev(y_k)) )
   * Toeplitz data flows from the bottom-right corner. We follow the data
   * by walking diagonally DOWN and RIGHT. Accumulate G[i+x] * H[j+x]
   * walking until we hit the bottom or right wall (max_x = min(n-1-i, m-1-j)).
   */

  int status = GR_SUCCESS;
  slong rank = gr_mat_ncols(G, ctx);

  if (gr_mat_nrows(A, ctx) != gr_mat_nrows(G, ctx) || gr_mat_ncols(A, ctx) != gr_mat_nrows(H, ctx)) { return GR_UNABLE; }

  gr_ptr sum_res = gr_heap_init(ctx); // final value for A[i,j]
  gr_ptr temp = gr_heap_init(ctx);    // g * h

  if (type == DISP_MINUS)
  {
    for (slong i = 0; i < gr_mat_nrows(G, ctx); i++)
    {
      for (slong j = 0; j < gr_mat_nrows(H, ctx); j++)
      {
        status |= gr_zero(sum_res, ctx);
        slong max_x = FLINT_MIN(gr_mat_nrows(G, ctx) - 1 - i, gr_mat_nrows(H, ctx) - 1 - j);
        for (slong k = 0; k < rank; k++)
        {
          for (slong x = 0; x <= max_x; x++)
          {
            status |= gr_mul(temp, gr_mat_entry_ptr(G, i + x, k, ctx), gr_mat_entry_ptr(H, j + x, k, ctx), ctx);
            status |= gr_add(sum_res, sum_res, temp, ctx);
          }
        }
        status |= gr_set(gr_mat_entry_ptr(A, i, j, ctx), sum_res, ctx);
      }
    }
  }
  else
  { // DISP_PLUS
    slong n = gr_mat_nrows(G, ctx);
    slong m = gr_mat_nrows(H, ctx);
    slong alpha = gr_mat_ncols(G, ctx);
    status = gr_mat_zero(A, ctx);
    if (status) return status;
    gr_poly_t pg, ph, p_tmp;
    for (slong k = 0; k < alpha; k++)
    {
      gr_poly_init(pg, ctx);
      gr_poly_fit_length(pg, n, ctx);
      for (slong r = 0; r < n; r++)
      {
        status |= gr_set(gr_poly_coeff_ptr(pg, r, ctx), gr_mat_entry_srcptr(G, r, k, ctx), ctx);
        if (status)
        {
          gr_poly_clear(pg, ctx);
          return status;
        }
      }
      _gr_poly_set_length(pg, n, ctx);
      for (slong j = 0; j < m; j++)
      {
        slong len = FLINT_MIN(j + 1, n);
        gr_poly_init(ph, ctx);
        gr_poly_fit_length(ph, len, ctx);
        for (slong r = 0; r < len; r++)
        {
          status |= gr_set(gr_poly_coeff_ptr(ph, r, ctx), gr_mat_entry_srcptr(H, j - r, k, ctx), ctx);
          if (status)
          {
            gr_poly_clear(pg, ctx);
            gr_poly_clear(ph, ctx);
            return status;
          }
        }
        _gr_poly_set_length(ph, len, ctx);
        gr_poly_init(p_tmp, ctx);
        status |= gr_poly_mul(p_tmp, pg, ph, ctx);
        if (status)
        {
          gr_poly_clear(pg, ctx);
          gr_poly_clear(ph, ctx);
          gr_poly_clear(p_tmp, ctx);
          return status;
        }
        if (gr_poly_length(p_tmp, ctx) > n) _gr_poly_set_length(p_tmp, n, ctx);
        for (slong i = 0; i < gr_poly_length(p_tmp, ctx); i++)
        {
          status |= gr_add(gr_mat_entry_ptr(A, i, j, ctx), gr_mat_entry_srcptr(A, i, j, ctx), gr_poly_coeff_srcptr(p_tmp, i, ctx), ctx);
          if (status)
          {
            gr_poly_clear(pg, ctx);
            gr_poly_clear(ph, ctx);
            gr_poly_clear(p_tmp, ctx);
            return status;
          }
        }
        gr_poly_clear(ph, ctx);
        gr_poly_clear(p_tmp, ctx);

        if (status)
        {
          gr_poly_clear(pg, ctx);
          return status;
        }
      }
      gr_poly_clear(pg, ctx);
    }
  }

  gr_heap_clear(sum_res, ctx);
  gr_heap_clear(temp, ctx);
  return status;
}