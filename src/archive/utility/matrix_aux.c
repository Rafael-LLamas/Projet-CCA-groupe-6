#include "matrix_aux.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"

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

// Application of the matrix Z on M, returnz to Res
int gr_mat_apply_Z(gr_mat_t Res, gr_mat_t M, gr_ctx_t ctx) {

  int error;

  slong n = gr_mat_nrows(M, ctx);

  slong m = gr_mat_ncols(M, ctx);

  for (slong j = 0; j < m; j++) {

    // On part du bas vers le haut pour ne pas écraser

    // les données avant qu'elles ne soient déplacées

    for (slong i = n - 1; i > 0; i--) {

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
int gr_mat_apply_Zt(gr_mat_t Res, gr_mat_t M, gr_ctx_t ctx) {

  int error;

  slong n = gr_mat_nrows(M, ctx);

  slong m = gr_mat_ncols(M, ctx);

  for (slong j = 0; j < m; j++) {

    // Pour Zt (vers le haut), on part du haut vers le bas

    for (slong i = 0; i < n - 1; i++) {

      error = gr_set(gr_mat_entry_ptr(Res, i, j, ctx), gr_mat_entry_srcptr(M, i + 1, j, ctx), ctx);

      if (error != 0) { return error; }
    }

    // La dernière ligne devient 0

    error = gr_zero(gr_mat_entry_ptr(Res, n - 1, j, ctx), ctx);

    if (error != 0) { return error; }
  }

  return GR_SUCCESS;
}