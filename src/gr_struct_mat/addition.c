#include <stdlib.h>
#include <time.h>

#include "flint/gr_mat.h"
#include "flint/gr_types.h"
#include "gr_struct_mat.h"

// Matrix addition of two structured matrix types
int gr_struct_mat_add(gr_struct_mat_t dest, gr_struct_mat_t mat, gr_ctx_t ctx)
{
  int status = GR_SUCCESS;

  if (gr_struct_mat_nrows(dest, ctx) != gr_struct_mat_nrows(mat, ctx) || gr_struct_mat_ncols(dest, ctx) != gr_struct_mat_ncols(mat, ctx))
    return GR_DOMAIN;

  // TODO: Future T_TOEPLITZ+HANKEL implementation
  if (dest->struct_t != mat->struct_t || dest->disp_t != mat->disp_t) return GR_DOMAIN;

  gr_mat_t temp_G, temp_H;
  gr_mat_init(temp_G, gr_mat_nrows(dest->G, ctx), gr_mat_ncols(dest->G, ctx) + gr_mat_ncols(mat->G, ctx), ctx);
  gr_mat_init(temp_H, gr_mat_nrows(dest->H, ctx), gr_mat_ncols(dest->H, ctx) + gr_mat_ncols(mat->H, ctx), ctx);

  status |= gr_mat_concat_horizontal(temp_G, dest->G, mat->G, ctx);
  status |= gr_mat_concat_horizontal(temp_H, dest->H, mat->H, ctx);

  if (status == GR_SUCCESS)
  {
    gr_mat_swap(dest->G, temp_G, ctx);
    gr_mat_swap(dest->H, temp_H, ctx);
  }

  gr_mat_clear(temp_G, ctx);
  gr_mat_clear(temp_H, ctx);
  return status;
}

// Matrix addition of a dense matrix to a structured matrix
int gr_struct_mat_add_dense(gr_struct_mat_t dest, gr_mat_t mat, structure_type_t struct_t, gr_ctx_t ctx)
{
  int status = GR_SUCCESS;
  gr_struct_mat_t s_mat;

  status |= gr_struct_mat_init_set(s_mat, mat, struct_t, dest->disp_t, ctx);
  status |= gr_struct_mat_add(dest, s_mat, ctx);
  gr_struct_mat_clear(s_mat, ctx);

  return status;
}
