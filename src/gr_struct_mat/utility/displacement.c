#include <stdio.h>
#include <stdlib.h>

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"
#include "gr_struct_mat.h"

// Takes the dense matrix A and returns the displaced matrix D
int _gr_mat_displacement(gr_mat_t D, gr_mat_t A, structure_type_t struct_t, disp_type_t type, gr_ctx_t ctx)
{
  /*
   * For Toeplitz matrices:
   * DISP_MINUS -> D = A - Z^T A Z   (A[i,j] - A[i+1, j+1])
   * DISP_PLUS  -> D = A - Z A Z^T   (A[i,j] - A[i-1, j-1])
   * * For Hankel matrices:
   * DISP_MINUS -> D = A - Z^T A Z^T (A[i,j] - A[i+1, j-1])
   * DISP_PLUS  -> D = A - Z A Z     (A[i,j] - A[i-1, j+1])
   */

  int status = GR_SUCCESS;
  slong n = gr_mat_nrows(A, ctx);
  slong m = gr_mat_ncols(A, ctx);

  if (gr_mat_nrows(D, ctx) != n || gr_mat_ncols(D, ctx) != m) return GR_UNABLE;
  if (n == 0 || m == 0) return GR_SUCCESS; // do nothing and success if 0x0 (we discussed this for too long)

  gr_ptr ptr_cur, ptr_op, ptr_dest;

  if (struct_t == T_TOEPLITZ)
  {
    if (type == DISP_MINUS)
    {
      for (slong i = 0; i < n - 1; i++)
        for (slong j = 0; j < m - 1; j++)
        {
          ptr_cur = gr_mat_entry_ptr(A, i, j, ctx);
          ptr_op = gr_mat_entry_ptr(A, i + 1, j + 1, ctx);
          ptr_dest = gr_mat_entry_ptr(D, i, j, ctx);
          status |= gr_sub(ptr_dest, ptr_cur, ptr_op, ctx);
        }
      // Copy last row and last column
      for (slong i = 0; i < n; i++) status |= gr_set(gr_mat_entry_ptr(D, i, m - 1, ctx), gr_mat_entry_ptr(A, i, m - 1, ctx), ctx);
      for (slong j = 0; j < m - 1; j++) status |= gr_set(gr_mat_entry_ptr(D, n - 1, j, ctx), gr_mat_entry_ptr(A, n - 1, j, ctx), ctx);
    }
    else // DISP_PLUS
    {
      for (slong i = 1; i < n; i++)
        for (slong j = 1; j < m; j++)
        {
          ptr_cur = gr_mat_entry_ptr(A, i, j, ctx);
          ptr_op = gr_mat_entry_ptr(A, i - 1, j - 1, ctx);
          ptr_dest = gr_mat_entry_ptr(D, i, j, ctx);
          status |= gr_sub(ptr_dest, ptr_cur, ptr_op, ctx);
        }
      // Copy last row and last column
      for (slong i = 0; i < n; i++) status |= gr_set(gr_mat_entry_ptr(D, i, 0, ctx), gr_mat_entry_ptr(A, i, 0, ctx), ctx);
      for (slong j = 1; j < m; j++) status |= gr_set(gr_mat_entry_ptr(D, 0, j, ctx), gr_mat_entry_ptr(A, 0, j, ctx), ctx);
    }
  }
  else if (struct_t == T_HANKEL)
  {
    if (type == DISP_MINUS)
    {
      for (slong i = 0; i < n - 1; i++)
      {
        for (slong j = 1; j < m; j++)
        {
          ptr_cur = gr_mat_entry_ptr(A, i, j, ctx);
          ptr_op = gr_mat_entry_ptr(A, i + 1, j - 1, ctx);
          ptr_dest = gr_mat_entry_ptr(D, i, j, ctx);
          status |= gr_sub(ptr_dest, ptr_cur, ptr_op, ctx);
        }
      }
      // Copy last row and first column
      for (slong j = 0; j < m; j++) status |= gr_set(gr_mat_entry_ptr(D, n - 1, j, ctx), gr_mat_entry_ptr(A, n - 1, j, ctx), ctx);
      for (slong i = 0; i < n - 1; i++) status |= gr_set(gr_mat_entry_ptr(D, i, 0, ctx), gr_mat_entry_ptr(A, i, 0, ctx), ctx);
    }
    else // DISP_PLUS
    {
      for (slong i = 1; i < n; i++)
      {
        for (slong j = 0; j < m - 1; j++)
        {
          ptr_cur = gr_mat_entry_ptr(A, i, j, ctx);
          ptr_op = gr_mat_entry_ptr(A, i - 1, j + 1, ctx);
          ptr_dest = gr_mat_entry_ptr(D, i, j, ctx);
          status |= gr_sub(ptr_dest, ptr_cur, ptr_op, ctx);
        }
      }
      // Copy first row and last column
      for (slong j = 0; j < m; j++) status |= gr_set(gr_mat_entry_ptr(D, 0, j, ctx), gr_mat_entry_ptr(A, 0, j, ctx), ctx);
      for (slong i = 1; i < n; i++) status |= gr_set(gr_mat_entry_ptr(D, i, m - 1, ctx), gr_mat_entry_ptr(A, i, m - 1, ctx), ctx);
    }
  }
  else
  {
    // Fail safely if structure type is T_UNSURE or unknown
    return GR_DOMAIN;
  }

  return status;
}