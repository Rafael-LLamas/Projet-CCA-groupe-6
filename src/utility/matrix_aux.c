#include "matrix_aux.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "random_toeplitz.h"

#include <stdlib.h>
#include <time.h>

/*
mets dans une matrices déja initialiser (voir pour init la matrice direct dans la fonction) des éléments random (<
limit) au position (i,j)
*/
int gr_mat_random_manual(gr_mat_t D, flint_rand_t state, ulong limit, gr_ctx_t ctx) {
  int error;
  for (int i = 0; i < gr_mat_nrows(D, ctx); i++) {
    for (int j = 0; j < gr_mat_ncols(D, ctx); j++) {
      error = gr_set_ui(gr_mat_entry_ptr(D, i, j, ctx), n_randint(state, limit), ctx);
      if (error != 0) return error;
    }
  }

  return GR_SUCCESS;
}

/*
La matrice LU doit avoir U et L avec la partie sup qui doit etre U sur M ligne et des 0 apres, L est toute la partie
triangulaire basse
*/
int gr_mat_lu_detach(gr_mat_t L, gr_mat_t U, gr_mat_t LU, gr_ctx_t ctx) {
  int error;
  for (slong i = 0; i < gr_mat_nrows(LU, ctx); i++) {
    error = gr_set_ui(gr_mat_entry_ptr(L, i, i, ctx), 1, ctx); // set 1 to diag.
    if (error != 0) return error;
    for (slong j = 0; j < gr_mat_ncols(LU, ctx); j++) { // distribute vals
      if (j >= i)                                       // diag and upper -> U
        error = gr_set(gr_mat_entry_ptr(U, i, j, ctx), gr_mat_entry_srcptr(LU, i, j, ctx), ctx);
      else // lower -> L
        error = gr_set(gr_mat_entry_ptr(L, i, j, ctx), gr_mat_entry_srcptr(LU, i, j, ctx), ctx);
      if (error != 0) return error;
    }
  }

  // for (int i = 0; i < gr_mat_ncols(LU, ctx); i++) {
  //   for (int j = i; j < gr_mat_ncols(LU, ctx); j++) {
  //     error = gr_set(gr_mat_entry_ptr(U, i, j, ctx), gr_mat_entry_srcptr(LU, i, j, ctx), ctx);
  //     if (error != 0) return error;
  //   }
  // }
  // for (int i = 0; i < gr_mat_nrows(LU, ctx); i++) {
  //   for (int j = 0; j < gr_mat_ncols(LU, ctx); j++) {
  //     if (i == j) {
  //       error = gr_set_ui(gr_mat_entry_ptr(L, i, j, ctx), 1, ctx);
  //       if (error != 0) return error;
  //     }
  //   }
  // }
  // for (int i = 1; i < gr_mat_nrows(LU, ctx); i++) {
  //   for (int j = 0; j < i; j++) {
  //     error = gr_set(gr_mat_entry_ptr(L, i, j, ctx), gr_mat_entry_srcptr(LU, i, j, ctx), ctx);
  //     if (error != 0) return error;
  //   }
  // }
  return GR_SUCCESS;
}