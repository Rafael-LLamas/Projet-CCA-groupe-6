#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"

#include <stdlib.h>
#include <time.h>

// Application of the matrix Z on M, returnz to Res
int _gr_mat_apply_Z(gr_mat_t Res, gr_mat_t M, gr_ctx_t ctx)
{

  int error;
  slong n = gr_mat_nrows(M, ctx);
  slong m = gr_mat_ncols(M, ctx);

  for (slong j = 0; j < m; j++)
  {
    // On part du bas vers le haut pour ne pas écraser
    // les données avant qu'elles ne soient déplacées
    for (slong i = n - 1; i > 0; i--)
    {
      error = gr_set(gr_mat_entry_ptr(Res, i, j, ctx), gr_mat_entry_srcptr(M, i - 1, j, ctx), ctx);
      if (error != 0) { return error; }
    }

    // La ligne 0 devient 0
    error = gr_zero(gr_mat_entry_ptr(Res, 0, j, ctx), ctx);
    if (error != 0) { return error; }
  }

  return GR_SUCCESS;
}

// Application of the transposed matrix Z on M, returnz to Res
int _gr_mat_apply_Zt(gr_mat_t Res, gr_mat_t M, gr_ctx_t ctx)
{
  int error;
  slong n = gr_mat_nrows(M, ctx);
  slong m = gr_mat_ncols(M, ctx);

  for (slong j = 0; j < m; j++)
  {
    // Pour Zt (vers le haut), on part du haut vers le bas
    for (slong i = 0; i < n - 1; i++)
    {
      error = gr_set(gr_mat_entry_ptr(Res, i, j, ctx), gr_mat_entry_srcptr(M, i + 1, j, ctx), ctx);
      if (error != 0) { return error; }
    }

    // La dernière ligne devient 0
    error = gr_zero(gr_mat_entry_ptr(Res, n - 1, j, ctx), ctx);
    if (error != 0) { return error; }
  }

  return GR_SUCCESS;
}