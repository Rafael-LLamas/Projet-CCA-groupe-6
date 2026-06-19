#include "../../gr_struct_mat.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"

#include <stdlib.h>
#include <time.h>

// Guess the structure type of the matrix mat, return it on type based off the accuracy rate in [0,1]
int _gr_struct_mat_guess_struct_t(structure_type_t type, gr_mat_t mat, slong acc_rate, gr_ctx_t ctx)
{
  int res = GR_SUCCESS;

  if (acc_rate > 1) acc_rate = 1;
  if (acc_rate < 0)
  {
    res = T_TOEPLITZ;
    return res;
  }

  slong nb_test = gr_mat_nrows(mat, ctx) * gr_mat_ncols(mat, ctx) * acc_rate;
  slong c_toeplitz = 0;
  slong c_hankel = 0;

  // https://arxiv.org/pdf/math-ph/0609050
  // Haar measure
  // Might wanna do something smart here
}