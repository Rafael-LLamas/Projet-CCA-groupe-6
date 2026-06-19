#include <stdio.h>
#include <stdlib.h>

#include "../../gr_struct_mat.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"


// Takes the dense matrix A and returns the displaced matrix D
int _gr_mat_displacement(gr_mat_t D, gr_mat_t A, disp_type_t type, gr_ctx_t ctx)
{

  /*
   * Computational approach to the previously seen displacement matrix generation.
   *
   * This function is an pointer arithmetic implementation of the previously
   * implemented gr_mat_displacement_square_safe(), we have 2 main options
   * for generation of two different displacement matrices mentioned as Phi-
   * and Phi+ on materials.
   *
   * For displacment operators:
   * Two types can be passed down to the function with the parameter type,
   * the type definitions are located in header file of this file.
   *
   *          DISP_MINUS  -> Phi-
   *          DISP_PLUS   -> Phi+ = ∇A
   *
   * Both will respectfull give the resulting matrix of A−Z^TAZ or A−ZAZ^T.
   */

  int status = GR_SUCCESS;
  slong n = gr_mat_nrows(A, ctx);
  slong m = gr_mat_ncols(A, ctx);

  if (gr_mat_nrows(D, ctx) != n || gr_mat_ncols(D, ctx) != m) return GR_UNABLE;
  if (n == 0 || m == 0) return GR_SUCCESS; // nothing to compute
  gr_ptr ptr_cur, ptr_op, ptr_dest;

  if (type == DISP_MINUS)
  {

    // inner part of Phi-
    for (slong i = 0; i < n - 1; i++)
    {
      for (slong j = 0; j < m - 1; j++)
      {
        ptr_cur = gr_mat_entry_ptr(A, i, j, ctx);
        ptr_op = gr_mat_entry_ptr(A, i + 1, j + 1, ctx);
        ptr_dest = gr_mat_entry_ptr(D, i, j, ctx);
        status |= gr_sub(ptr_dest, ptr_cur, ptr_op, ctx); // D[i,j] = A[i,j] - A[i+1, j+1]
      }
    }

    // copy the rest (last row & column)
    for (slong i = 0; i < n; i++) status |= gr_set(gr_mat_entry_ptr(D, i, m - 1, ctx), gr_mat_entry_ptr(A, i, m - 1, ctx), ctx);
    for (slong j = 1; j < m; j++) status |= gr_set(gr_mat_entry_ptr(D, n - 1, j, ctx), gr_mat_entry_ptr(A, n - 1, j, ctx), ctx);
  }
  else
  { // DISP_PLUS

    // inner part of Phi+
    for (slong i = n - 1; i > 0; i--)
    {
      for (slong j = 1; j < m; j++)
      {
        ptr_cur = gr_mat_entry_ptr(A, i, j, ctx);
        ptr_op = gr_mat_entry_ptr(A, i - 1, j - 1, ctx);
        ptr_dest = gr_mat_entry_ptr(D, i, j, ctx);
        status |= gr_sub(ptr_dest, ptr_cur, ptr_op, ctx); // D[i,j] = A[i,j] - A[i-1, j-1]
      }
    }

    // copy the rest (first row & columns)
    for (slong i = 0; i < n; i++) status |= gr_set(gr_mat_entry_ptr(D, i, 0, ctx), gr_mat_entry_ptr(A, i, 0, ctx), ctx);
    for (slong j = 1; j < m; j++) status |= gr_set(gr_mat_entry_ptr(D, 0, j, ctx), gr_mat_entry_ptr(A, 0, j, ctx), ctx);
  }

  return status;
}