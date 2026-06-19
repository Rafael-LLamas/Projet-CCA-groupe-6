#include "../gr_struct_mat.h"
#include "flint/flint.h"
#include "flint/gr_types.h"

#include "utility/utility.h"

// Initialize the structured matrix and allocate memory for G and H
int gr_struct_mat_init_set(gr_struct_mat_t struct_mat, gr_mat_t mat, structure_type_t mat_struct, disp_type_t disp, gr_ctx_t ctx)
{
  int res = GR_SUCCESS;

  

  
  res |= _gr_mat_G_H(gr_mat_struct * G, gr_mat_struct * H, gr_mat_struct * A, disp_type_t type, gr_ctx_struct * ctx) struct_mat->struct_t = mat;
  mat->disp_t = disp;
  gr_mat_init(mat->G, n, m, ctx);
  gr_mat_init(mat->H, n, m, ctx);
}

// Free the memory of G and H
void gr_struct_mat_clear(gr_struct_mat_t mat, gr_ctx_t ctx)
{
  gr_mat_clear(mat->G, ctx);
  gr_mat_clear(mat->H, ctx);
}

// Return the number of rows
slong gr_struct_mat_nrows(gr_struct_mat_srcptr mat, gr_ctx_t ctx) { return gr_mat_nrows(mat->G, ctx); }

// Return the rank of the matrix (from generators)
slong gr_struct_mat_rank(gr_struct_mat_srcptr mat, gr_ctx_t ctx) { return gr_mat_ncols(mat->G, ctx); }

// Print the structured matrix with additional information
void gr_struct_mat_print(gr_struct_mat_srcptr mat, gr_ctx_t ctx)
{
  flint_printf("Structured Matrix: %s, Displacement: %s\n", mat->struct_t == TOEPLITZ ? "Toeplitz" : "Hankel",
               mat->disp_t == DISP_PLUS ? "Phi_Plus" : "Phi_Minus");
  flint_printf("--- Generator G ---\n");
  gr_mat_print(mat->G, ctx);
  flint_printf("\n--- Generator H ---\n");
  gr_mat_print(mat->H, ctx);
  flint_printf("\n");
}

int gr_struct_mat_