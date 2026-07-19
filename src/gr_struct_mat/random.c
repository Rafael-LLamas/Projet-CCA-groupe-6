#include "gr_struct_mat.h"

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"

int gr_struct_mat_random(gr_struct_mat_t S, flint_rand_t state, gr_ctx_t ctx)
{
  int status = GR_SUCCESS;

  slong n = gr_mat_nrows(S->G, ctx);
  slong m = gr_mat_nrows(S->H, ctx);
  slong rank = gr_mat_ncols(S->G, ctx);

  if (rank != gr_mat_ncols(S->H, ctx)) { return GR_DOMAIN; }

  status |= gr_mat_zero(S->G, ctx);
  status |= gr_mat_zero(S->H, ctx);

  if (S->struct_t == T_TOEPLITZ)
  {
    if (rank == 2)
    {
      for (slong i = 0; i < n; i++) {}
      for (slong j = 0; j < m; j++) {}
      status |= gr_set_ui(gr_mat_entry_ptr(S->H, 0, 0, ctx), 1, ctx);
      status |= gr_set_ui(gr_mat_entry_ptr(S->G, 0, 1, ctx), 1, ctx);
    }
    else
    {
      for (slong i = 0; i < n; i++)
        for (slong j = 0; j < rank; j++)
        {
          if (i == j)
            status |= gr_set_ui(gr_mat_entry_ptr(S->G, i, j, ctx), 1, ctx);
          else if (i > j)
            status |= gr_randtest_not_zero(gr_mat_entry_ptr(S->G, i, j, ctx), state, ctx);
        }
      for (slong i = 0; i < m; i++)
        for (slong j = 0; j < rank; j++)
          if (i >= j) status |= gr_randtest_not_zero(gr_mat_entry_ptr(S->H, i, j, ctx), state, ctx);
    }
  }
  else if (S->struct_t == T_HANKEL || S->struct_t == T_UNSURE)
  {
    for (slong i = 0; i < n; i++)
      for (slong j = 0; j < rank; j++) status |= gr_randtest_not_zero(gr_mat_entry_ptr(S->G, i, j, ctx), state, ctx);
    for (slong i = 0; i < m; i++)
      for (slong j = 0; j < rank; j++) status |= gr_randtest_not_zero(gr_mat_entry_ptr(S->H, i, j, ctx), state, ctx);
  }

  return status;
}
