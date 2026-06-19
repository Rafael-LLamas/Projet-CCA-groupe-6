#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"

#include <stdlib.h>
#include <time.h>

/*
La matrice LU doit avoir U et L avec la partie sup qui doit etre U sur M ligne et des 0 apres, L est toute la partie
triangulaire basse
*/
int _gr_mat_lu_detach(gr_mat_t L, gr_mat_t U, gr_mat_t LU, gr_ctx_t ctx)
{
  int error;
  for (slong i = 0; i < gr_mat_nrows(LU, ctx); i++)
  {
    error = gr_set_ui(gr_mat_entry_ptr(L, i, i, ctx), 1, ctx); // set 1 to diag.
    if (error != 0) return error;
    for (slong j = 0; j < gr_mat_ncols(LU, ctx); j++)
    {             // distribute vals
      if (j >= i) // diag and upper -> U
        error = gr_set(gr_mat_entry_ptr(U, i, j, ctx), gr_mat_entry_srcptr(LU, i, j, ctx), ctx);
      else // lower -> L
        error = gr_set(gr_mat_entry_ptr(L, i, j, ctx), gr_mat_entry_srcptr(LU, i, j, ctx), ctx);
      if (error != 0) return error;
    }
  }
  return GR_SUCCESS;
}