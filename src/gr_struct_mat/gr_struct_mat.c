#include "gr_struct_mat.h"
#include "flint/flint.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"

#include "utility/utility.h"

// Initialize the structured matrix shell
void gr_struct_mat_init(gr_struct_mat_t struct_mat, slong rows, slong cols, structure_type_t mat_struct, disp_type_t disp_t, gr_ctx_t ctx)
{
  struct_mat->struct_t = mat_struct;
  struct_mat->disp_t = disp_t;
  gr_mat_init(struct_mat->G, rows, 0, ctx);
  gr_mat_init(struct_mat->H, cols, 0, ctx);
}

// Initialize the structured matrix from an existing dense matrix
int gr_struct_mat_init_set(gr_struct_mat_t struct_mat, gr_mat_t mat, structure_type_t mat_struct, disp_type_t disp_t, gr_ctx_t ctx)
{
  int status = GR_SUCCESS;
  if (mat_struct == T_UNSURE) { status |= _gr_struct_mat_guess_struct_t(&mat_struct, mat, 1, ctx); }
  struct_mat->disp_t = disp_t;
  struct_mat->struct_t = mat_struct;
  status |= _gr_mat_G_H(struct_mat->G, struct_mat->H, mat, disp_t, ctx);
  return status;
}

// Free the memory of G and H
void gr_struct_mat_clear(gr_struct_mat_t mat, gr_ctx_t ctx)
{
  gr_mat_clear(mat->G, ctx);
  gr_mat_clear(mat->H, ctx);
}

// Return the number of rows
slong gr_struct_mat_nrows(gr_struct_mat_t mat, gr_ctx_t ctx) { return gr_mat_nrows(mat->G, ctx); }

// Return the number of cols
slong gr_struct_mat_ncols(gr_struct_mat_t mat, gr_ctx_t ctx) { return gr_mat_nrows(mat->H, ctx); }

// Return the rank of the matrix (from generators)
slong gr_struct_mat_rank(gr_struct_mat_t mat, gr_ctx_t ctx)
{
  return (gr_mat_ncols(mat->G, ctx) < gr_mat_ncols(mat->H, ctx)) ? gr_mat_ncols(mat->G, ctx) : gr_mat_ncols(mat->H, ctx);
}

// Print the structured matrix with additional information
void gr_struct_mat_print(gr_struct_mat_t mat, gr_ctx_t ctx)
{
  flint_printf("Structured Matrix: %s, Displacement: %s\n", mat->struct_t == T_TOEPLITZ ? "Toeplitz" : "Hankel",
               mat->disp_t == DISP_PLUS ? "Phi_Plus" : "Phi_Minus");
  flint_printf("--- Generator G ---\n");
  gr_mat_print(mat->G, ctx);
  flint_printf("\n--- Generator H ---\n");
  gr_mat_print(mat->H, ctx);
  flint_printf("\n");
}

int gr_struct_mat_reconstruct(gr_mat_t dense_mat, gr_struct_mat_t mat, gr_ctx_t ctx)
{
  int status = GR_SUCCESS;

  if (gr_mat_nrows(dense_mat, ctx) != gr_struct_mat_nrows(mat, ctx) || gr_mat_ncols(dense_mat, ctx) != gr_struct_mat_ncols(mat, ctx))
  {
    gr_mat_clear(dense_mat, ctx);
    gr_mat_init(dense_mat, gr_struct_mat_nrows(mat, ctx), gr_struct_mat_ncols(mat, ctx), ctx);
  }

  status |= _gr_mat_reconstruct(dense_mat, mat->G, mat->H, mat->disp_t, ctx);

  return status;
}