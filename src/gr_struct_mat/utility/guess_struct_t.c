#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "gr_struct_mat.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"

int _gr_struct_mat_guess_struct_t(structure_type_t *type, gr_mat_t mat, float acc_rate, gr_ctx_t ctx)
{
  int res = GR_SUCCESS;

  if (type == NULL) return GR_DOMAIN;

  if (acc_rate > 1) acc_rate = 1.0f;
  if (acc_rate < 0)
  {
    *type = T_TOEPLITZ;
    return res;
  }

  slong rows = gr_mat_nrows(mat, ctx);
  slong cols = gr_mat_ncols(mat, ctx);

  // if matrix is 1x1 then assume Toeplitz but no computational advantage possible
  if (rows <= 1 || cols <= 1)
  {
    *type = T_TOEPLITZ;
    return res;
  }

  // A better approach possible on checking only req_check elements rather than O(n^2)
  // Careful on selecting "uniformly" random elemnets
  // https://arxiv.org/pdf/math-ph/0609050
  // Haar measure?

  slong nb_check = (rows - 1) * (cols - 1);
  slong req_check = (slong)round((double)nb_check * acc_rate);

  slong c_toeplitz = 0;

  for (slong i = 0; i < rows - 1; i++)
    for (slong j = 0; j < cols - 1; j++)
    {
      gr_ptr toep_curr = gr_mat_entry_ptr(mat, i, j, ctx);
      gr_ptr toep_next = gr_mat_entry_ptr(mat, i + 1, j + 1, ctx);
      if (gr_equal(toep_curr, toep_next, ctx) == T_TRUE) { c_toeplitz++; }
    }

  if (c_toeplitz >= req_check)
    *type = T_TOEPLITZ;
  else
    *type = T_HANKEL;

  return res;
}