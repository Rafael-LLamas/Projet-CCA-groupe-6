#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "matrix_aux.h"
#include "random_toeplitz.h"

#include <stdlib.h>
#include <time.h>

int gr_multiplication_toeplitz(gr_mat_t C, gr_mat_t A, gr_mat_t B, gr_ctx_t ctx) {
  /*
  Si j'ai bien tout compris, on doit d'abord faire sur C le produit scalaire complet pour la ligne 0 et la colonne 0.
  Puis je prends en compte que la 1ere colonne et ligne  de A et B.
  Avec ca je dois suivre la formule : C i,j = Ci-1,j-1 + (Ai,n-1 x Bn-1,j) - (Ai-1,0 x B0,j-1), sachant que les matrices
  sont Toeplitz donc je peux simplifier A et B par leurs vecteur et juste prendre la bonne position
  */
  int error;
  for (int i = 0; i < gr_mat_nrows(C, ctx); i++) {
    for (int j = 0; j < gr_mat_ncols(A, ctx); j++) {
      error = gr_mul(gr_mat_entry_ptr(C, i, 0, ctx), gr_mat_entry_srcptr(A, i, j, ctx),
                     gr_mat_entry_srcptr(B, j, 0, ctx), ctx);
    }
  }
  for (int j = 1; j < gr_mat_ncols(C, ctx); j++) {
    for (int i = 0; i < gr_mat_ncols(A, ctx); i++) {
      error = gr_mul(gr_mat_entry_ptr(C, 0, j, ctx), gr_mat_entry_srcptr(A, 0, i, ctx),
                     gr_mat_entry_srcptr(B, i, j, ctx), ctx);
    }
  }
  gr_ptr temp;
  GR_TMP_INIT(temp, ctx);
  for (int i = 1; i < gr_mat_nrows(C, ctx); i++) {
    for (int j = 1; j < gr_mat_ncols(C, ctx); j++) {
      error = gr_mul(gr_mat_entry_ptr(C, i, j, ctx), gr_mat_entry_srcptr(A, i, gr_mat_ncols(A, ctx) - 1, ctx),
                     gr_mat_entry_srcptr(B, gr_mat_ncols(A, ctx) - 1, j, ctx), ctx);
      error = gr_add(gr_mat_entry_ptr(C, i, j, ctx), gr_mat_entry_srcptr(C, i - 1, j - 1, ctx),
                     gr_mat_entry_ptr(C, i, j, ctx), ctx);
      gr_mul(temp, gr_mat_entry_srcptr(A, i - 1, 0, ctx), gr_mat_entry_srcptr(B, 0, j - 1, ctx), ctx);
      error = gr_sub(gr_mat_entry_ptr(C, i, j, ctx), gr_mat_entry_srcptr(C, i - 1, j - 1, ctx), (gr_srcptr)temp, ctx);
    }
  }
  return 0;
}