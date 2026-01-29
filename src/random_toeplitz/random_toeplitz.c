#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include <stdlib.h>

void random_quasi_toeplitz(gr_mat_t A, int n, int m, flint_rand_t state, gr_ctx_t ctx) {
  /*
  Mon sommeil va en patir ......
  Donc apres avoir compris mon problème ici je génère des quasi toeplitz en suivant la logique que mets matrices auront
  un rang faible. Pour ca je génère 2 matrice, L et U puis je l'ai multiplie entre elle pour avoir A.
  */
  gr_mat_t L, U, colL, rowU;
  int error;
  gr_mat_init(L, n, m, ctx);
  gr_mat_init(U, n, m, ctx);
  gr_mat_init(colL, n, 1, ctx);
  gr_mat_init(rowU, 1, m, ctx);
  gr_mat_init(A, n, m, ctx);
  error = gr_mat_zero(L, ctx);
  error = gr_mat_zero(U, ctx);

  for (int i = 0; i < n; i++) error = gr_randtest_not_zero(gr_mat_entry_ptr(colL, i, 0, ctx), state, ctx);
  for (int j = 0; j < m; j++) error = gr_randtest_not_zero(gr_mat_entry_ptr(rowU, 0, j, ctx), state, ctx);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      if (i >= j) { error = gr_set(gr_mat_entry_ptr(L, i, j, ctx), gr_mat_entry_srcptr(colL, i - j, 0, ctx), ctx); }
      if (j >= i) { error = gr_set(gr_mat_entry_ptr(U, i, j, ctx), gr_mat_entry_srcptr(rowU, 0, j - i, ctx), ctx); }
    }
  }
  error = gr_mat_mul(A, L, U, ctx);
  gr_mat_clear(L, ctx);
  gr_mat_clear(U, ctx);
  gr_mat_clear(colL, ctx);
  gr_mat_clear(rowU, ctx);
}

void random_toeplitz(gr_mat_t A, int n, int m, flint_rand_t state, gr_ctx_t ctx) {
  /*Ici je génère de façon non-chalante une matrice toeplitz en prenant 2 vecteur aléatoire que j'utilise pour def la
   * matrice*/
  gr_mat_t col, row;
  int error;
  gr_mat_init(col, n, 1, ctx);
  gr_mat_init(row, 1, m, ctx);

  error = gr_mat_randtest(col, state, ctx);
  if (error != 0) { exit(error); }
  error = gr_mat_randtest(row, state, ctx);
  if (error != 0) { exit(error); }
  error = gr_set(gr_mat_entry_ptr(row, 0, 0, ctx), gr_mat_entry_ptr(col, 0, 0, ctx), ctx);
  if (error != 0) { exit(error); }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      if (i >= j) {
        error = gr_set(gr_mat_entry_ptr(A, i, j, ctx), gr_mat_entry_srcptr(col, i - j, 0, ctx), ctx);
        if (error != 0) { exit(error); }
      } else {
        error = gr_set(gr_mat_entry_ptr(A, i, j, ctx), gr_mat_entry_srcptr(row, 0, j - i, ctx), ctx);
        if (error != 0) { exit(error); }
      }
    }
  }

  gr_mat_clear(col, ctx);
  gr_mat_clear(row, ctx);
}

int test_random_toeplitz() {
  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, 29);
  gr_mat_t ran;
  flint_rand_t state;
  flint_rand_init(state);
  gr_mat_init(ran, 5, 5, ctx);
  random_quasi_toeplitz(ran, 5, 5, state, ctx);
  flint_printf("----------------------------------------------\nResultat quasi = ");
  gr_mat_print(ran, ctx);
  flint_printf("\n");

  gr_mat_clear(ran, ctx);
  gr_ctx_clear(ctx);
  return GR_SUCCESS;
}