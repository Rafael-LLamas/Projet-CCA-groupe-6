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

// TODO - might integrate with cmake tests later
int test_matrix_aux() {
  flint_printf("*----------* Matrix Auxilary Test *----------*\n");
  gr_ctx_t ctx;
  {
    flint_printf("\n:-------: LU Detatch on random Matrix :-------:\n");
    gr_mat_t A, LU, L, U, B, C;
    slong rank, *P;
    flint_rand_t state;
    flint_rand_init(state);
    gr_ctx_init_nmod(ctx, 1009);
    P = flint_malloc(sizeof(slong) * 10);
    gr_mat_init(A, 5, 10, ctx);
    gr_mat_init(C, 10, 10, ctx);
    int error = gr_mat_init_set(B, A, ctx);
    if (error != 0) {
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      flint_rand_clear(state);
      flint_free(P);
      gr_ctx_clear(ctx);
      return error;
    }
    gr_mat_init(LU, 5, 10, ctx);
    gr_mat_init(L, 5, 10, ctx);
    gr_mat_init(U, 10, 10, ctx);
    error = rand_quasi_toeplitz(A, 5, 10, 0, ctx);
    if (error != 0) {
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_mat_clear(LU, ctx);
      gr_mat_clear(L, ctx);
      gr_mat_clear(U, ctx);
      flint_rand_clear(state);
      flint_free(P);
      gr_ctx_clear(ctx);
      return error;
    }
    flint_printf("A = \n");
    gr_mat_print(A, ctx);
    flint_printf("\n**************************************************\n");
    flint_printf("LU = \n");
    error = gr_mat_lu(&rank, P, LU, A, 0, ctx);
    if (error != 0) {
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_mat_clear(LU, ctx);
      gr_mat_clear(L, ctx);
      gr_mat_clear(U, ctx);
      flint_rand_clear(state);
      flint_free(P);
      gr_ctx_clear(ctx);
      return error;
    }
    gr_mat_print(LU, ctx);
    flint_printf("\n**************************************************\n");
    error = gr_mat_lu_detach(L, U, LU, ctx);
    if (error != 0) {
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_mat_clear(LU, ctx);
      gr_mat_clear(L, ctx);
      gr_mat_clear(U, ctx);
      flint_rand_clear(state);
      flint_free(P);
      gr_ctx_clear(ctx);
      return error;
    }
    flint_printf("U = \n");
    gr_mat_print(U, ctx);
    flint_printf("\n\n");
    flint_printf("L = \n");
    gr_mat_print(L, ctx);
    flint_printf("\n**************************************************\n");
    error = gr_mat_mul(B, L, U, ctx);
    if (error != 0) {
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_mat_clear(LU, ctx);
      gr_mat_clear(L, ctx);
      gr_mat_clear(U, ctx);
      flint_rand_clear(state);
      flint_free(P);
      gr_ctx_clear(ctx);
      return error;
    }
    flint_printf("B = \n");
    gr_mat_print(B, ctx);
    flint_printf("\n**************************************************\n");
    error = gr_mat_random_manual(C, state, 1009, ctx);
    if (error != 0) {
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_mat_clear(LU, ctx);
      gr_mat_clear(L, ctx);
      gr_mat_clear(U, ctx);
      flint_rand_clear(state);
      flint_free(P);
      gr_ctx_clear(ctx);
      return error;
    }
    flint_printf("C = \n");
    gr_mat_print(C, ctx);
    gr_mat_clear(A, ctx);
    gr_mat_clear(LU, ctx);
    gr_mat_clear(L, ctx);
    gr_mat_clear(U, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);
    flint_rand_clear(state);
    flint_free(P);
  }
  {
    flint_printf("\n\n:-------: LU Detatch on weird Matrix :-------:\n");
    gr_mat_t A, LU, L, U, B;
    slong rank, *P;
    gr_ctx_init_nmod(ctx, GNMOD);
    P = flint_malloc(sizeof(slong) * 10);
    gr_mat_init(A, 5, 4, ctx);
    gr_mat_init(LU, 5, 4, ctx);
    gr_mat_init(L, 5, 5, ctx);
    gr_mat_init(U, 5, 4, ctx);
    gr_mat_init(B, 5, 4, ctx);
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 0, ctx), 1, ctx)); // Row 0: 1, 2, 1, 2, 1
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 1, ctx), 2, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 2, ctx), 1, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 3, ctx), 2, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 4, ctx), 1, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 0, ctx), 2, ctx)); // Row 1: 2, 4, 3, 6, 1
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 1, ctx), 4, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 2, ctx), 3, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 3, ctx), 6, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 4, ctx), 1, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 0, ctx), 3, ctx)); // Row 2: 3, 6, 3, 6, 5
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 1, ctx), 6, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 2, ctx), 3, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 3, ctx), 6, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 4, ctx), 5, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 3, 0, ctx), 4, ctx)); // Row 3: 4, 8, 4, 8, 6
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 3, 1, ctx), 8, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 3, 2, ctx), 4, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 3, 3, ctx), 8, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 3, 4, ctx), 6, ctx));
    // FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 4, 0, ctx), 5, ctx)); // Row 4:  5, 10, 5, 10, 5,
    // FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 4, 1, ctx), 10, ctx));
    // FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 4, 2, ctx), 5, ctx));
    // FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 4, 3, ctx), 10, ctx));
    // FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 4, 4, ctx), 5, ctx));
    flint_printf("A = \n");
    gr_mat_print(A, ctx);
    flint_printf("\n**************************************************\n");
    flint_printf("LU = \n");
    FLINT_CHECK(gr_mat_lu(&rank, P, LU, A, 0, ctx));
    gr_mat_print(LU, ctx);
    flint_printf("\n**************************************************\n");
    FLINT_CHECK(gr_mat_lu_detach(L, U, LU, ctx));
    flint_printf("U = \n");
    gr_mat_print(U, ctx);
    flint_printf("\n\n");
    flint_printf("L = \n");
    gr_mat_print(L, ctx);
    flint_printf("\n**************************************************\n");
    FLINT_CHECK(gr_mat_mul(B, L, U, ctx));
    flint_printf("B = \n");
    gr_mat_print(B, ctx);
    flint_printf("\n");
    // for (slong i = 0; i < gr_mat_nrows(B, ctx); i++) FLINT_CHECK(gr_mat_move_row(B, i, P[i], ctx)); // this fails
    flint_printf("B but (not yet) permutated = \n");
    gr_mat_print(B, ctx);
    flint_printf("\n");
    gr_mat_clear(A, ctx);
    gr_mat_clear(LU, ctx);
    gr_mat_clear(L, ctx);
    gr_mat_clear(U, ctx);
    gr_mat_clear(B, ctx);
    flint_free(P);
  }
  gr_ctx_clear(ctx);
  return GR_SUCCESS;
}