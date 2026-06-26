#include <stdio.h>
#include <stdlib.h>

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"
#include "gr_struct_mat.h"

int _gr_mat_reconstruct_toeplitz(gr_mat_t A, gr_mat_t G, gr_mat_t H, disp_type_t type, gr_ctx_t ctx)
{
  int status = GR_SUCCESS;
  slong n = gr_mat_nrows(G, ctx);
  slong m = gr_mat_nrows(H, ctx);
  slong rank = gr_mat_ncols(G, ctx);

  if (gr_mat_nrows(A, ctx) != n || gr_mat_ncols(A, ctx) != m) return GR_UNABLE;

  gr_ptr dot_prod = gr_heap_init(ctx);
  gr_ptr temp = gr_heap_init(ctx);

  if (type == DISP_PLUS) // flow is top-left to bottom-right:
    for (slong i = 0; i < n; i++)
      for (slong j = 0; j < m; j++)
      {
        status |= gr_zero(dot_prod, ctx);
        for (slong k = 0; k < rank; k++)
        {
          status |= gr_mul(temp, gr_mat_entry_ptr(G, i, k, ctx), gr_mat_entry_ptr(H, j, k, ctx), ctx);
          status |= gr_add(dot_prod, dot_prod, temp, ctx);
        }

        if (i > 0 && j > 0)
          status |= gr_add(gr_mat_entry_ptr(A, i, j, ctx), gr_mat_entry_ptr(A, i - 1, j - 1, ctx), dot_prod, ctx);
        else
          status |= gr_set(gr_mat_entry_ptr(A, i, j, ctx), dot_prod, ctx);
      }
  else // DISP_MINUS
  {
    // flow bottom-right to top-left:
    for (slong i = n - 1; i >= 0; i--)
      for (slong j = m - 1; j >= 0; j--)
      {
        status |= gr_zero(dot_prod, ctx);
        for (slong k = 0; k < rank; k++)
        {
          status |= gr_mul(temp, gr_mat_entry_ptr(G, i, k, ctx), gr_mat_entry_ptr(H, j, k, ctx), ctx);
          status |= gr_add(dot_prod, dot_prod, temp, ctx);
        }

        if (i < n - 1 && j < m - 1)
          status |= gr_add(gr_mat_entry_ptr(A, i, j, ctx), gr_mat_entry_ptr(A, i + 1, j + 1, ctx), dot_prod, ctx);
        else
          status |= gr_set(gr_mat_entry_ptr(A, i, j, ctx), dot_prod, ctx);
      }
  }

  gr_heap_clear(dot_prod, ctx);
  gr_heap_clear(temp, ctx);
  return status;
}


int _gr_mat_reconstruct_hankel(gr_mat_t A, gr_mat_t G, gr_mat_t H, disp_type_t type, gr_ctx_t ctx)
{
  int status = GR_SUCCESS;
  slong n = gr_mat_nrows(G, ctx);
  slong m = gr_mat_nrows(H, ctx);
  slong rank = gr_mat_ncols(G, ctx);

  if (gr_mat_nrows(A, ctx) != n || gr_mat_ncols(A, ctx) != m) { return GR_UNABLE; }

  gr_ptr dot_prod = gr_heap_init(ctx);
  gr_ptr temp = gr_heap_init(ctx);

  if (type == DISP_PLUS) // flow is top-right to bottom-left:
    for (slong i = 0; i < n; i++)
      for (slong j = 0; j < m; j++)
      {
        status |= gr_zero(dot_prod, ctx);
        for (slong k = 0; k < rank; k++)
        {
          status |= gr_mul(temp, gr_mat_entry_ptr(G, i, k, ctx), gr_mat_entry_ptr(H, j, k, ctx), ctx);
          status |= gr_add(dot_prod, dot_prod, temp, ctx);
        }

        if (i > 0 && j < m - 1)
          status |= gr_add(gr_mat_entry_ptr(A, i, j, ctx), gr_mat_entry_ptr(A, i - 1, j + 1, ctx), dot_prod, ctx);
        else
          status |= gr_set(gr_mat_entry_ptr(A, i, j, ctx), dot_prod, ctx);
      }
  else // DISP_MINUS
  {
    // flow is bottom-left to top-right:
    for (slong i = n - 1; i >= 0; i--)
      for (slong j = 0; j < m; j++)
      {
        status |= gr_zero(dot_prod, ctx);
        for (slong k = 0; k < rank; k++)
        {
          status |= gr_mul(temp, gr_mat_entry_ptr(G, i, k, ctx), gr_mat_entry_ptr(H, j, k, ctx), ctx);
          status |= gr_add(dot_prod, dot_prod, temp, ctx);
        }

        if (i < n - 1 && j > 0)
          status |= gr_add(gr_mat_entry_ptr(A, i, j, ctx), gr_mat_entry_ptr(A, i + 1, j - 1, ctx), dot_prod, ctx);
        else
          status |= gr_set(gr_mat_entry_ptr(A, i, j, ctx), dot_prod, ctx);
      }
  }

  gr_heap_clear(dot_prod, ctx);
  gr_heap_clear(temp, ctx);
  return status;
}