#include "matrix_aux.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "random_toeplitz.h"

#include <stdlib.h>
#include <time.h>

int gr_mat_random_manual(gr_mat_t D, gr_ctx_t ctx) { return gr_not_implemented(); }

/*
La matrice LU doit avoir U et L avec la partie sup qui doit etre U sur M ligne et des 0 apres, L est toute la partie
triangulaire basse
*/
int gr_mat_lu_detach(gr_mat_t L, gr_mat_t U, gr_mat_t LU, gr_ctx_t ctx) {
  for (int i = 0; i < gr_mat_ncols(LU, ctx); i++) {
    for (int j = i; j < gr_mat_ncols(LU, ctx); j++) {
      gr_set(gr_mat_entry_ptr(U, i, j, ctx), gr_mat_entry_srcptr(LU, i, j, ctx), ctx);
    }
  }
  for (int i = 0; i < gr_mat_nrows(LU, ctx); i++) {
    for (int j = 0; j < gr_mat_ncols(LU, ctx); j++) {
      if (i == j) { gr_mat_set_ui(gr_mat_entry_ptr(L, i, j, ctx), 1, ctx); }
    }
  }
  for (int i = 1; i < gr_mat_nrows(LU, ctx); i++) {
    for (int j = 0; j < i; j++) { gr_set(gr_mat_entry_ptr(L, i, j, ctx), gr_mat_entry_srcptr(LU, i, j, ctx), ctx); }
  }
  return gr_not_implemented();
}

// TODO - might integrate with cmake tests later
int test_matrix_aux() {
  flint_printf("*----------* Matrix Auxilary Test *----------*\n");
  gr_ctx_t ctx;
  gr_mat_t A, LU, L, U, B;
  slong rank, *P;
  gr_ctx_init_nmod(ctx, 1009);
  P = flint_malloc(sizeof(slong) * 10);
  gr_mat_init(A, 5, 10, ctx);
  gr_mat_init_set(B, A, ctx);
  gr_mat_init(LU, 5, 10, ctx);
  gr_mat_init(L, 5, 10, ctx);
  gr_mat_init(U, 10, 10, ctx);
  int error = rand_quasi_toeplitz(A, 5, 10, 0, ctx);
  flint_printf("A = \n");
  gr_mat_print(A, ctx);
  flint_printf("\n**************************************************\n");
  flint_printf("LU = \n");
  error = gr_mat_lu(&rank, P, LU, A, 0, ctx);
  gr_mat_print(LU, ctx);
  flint_printf("\n**************************************************\n");
  int res = gr_mat_lu_detach(L, U, LU, ctx);
  flint_printf("U = \n");
  gr_mat_print(U, ctx);
  flint_printf("\n\n");
  flint_printf("L = \n");
  gr_mat_print(L, ctx);
  flint_printf("\n**************************************************\n");
  gr_mat_mul(B, L, U, ctx);
  flint_printf("B = \n");
  gr_mat_print(B, ctx);
  flint_printf("\n");
  gr_mat_clear(A, ctx);
  gr_mat_clear(LU, ctx);
  return res;
}