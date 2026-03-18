#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include <stdlib.h>

int gr_mat_random_quasi_toepitz(gr_mat_t A, int n, int m, flint_rand_t state, gr_ctx_t ctx) {
  /*
  Mon sommeil va en patir ......
  Donc apres avoir compris mon problème ici je génère des quasi toeplitz en suivant la logique que mets matrices auront
  un rang faible. Pour ca je génère 2 matrice toeplitz, L et U puis je l'ai multiplie entre elle pour avoir A.
  */
  gr_mat_t L, U, colL, rowU;
  int error;
  gr_mat_init(L, n, m, ctx);
  gr_mat_init(U, m, m, ctx);
  gr_mat_init(colL, n, 1, ctx);
  gr_mat_init(rowU, 1, m, ctx);
  error = gr_mat_zero(L, ctx);
  if (error != 0) { return error; }
  error = gr_mat_zero(U, ctx);
  if (error != 0) { return error; }

  for (int i = 0; i < n; i++) {
    error = gr_randtest_not_zero(gr_mat_entry_ptr(colL, i, 0, ctx), state, ctx);
    if (error != 0) { return error; }
  }
  for (int j = 0; j < m; j++) {
    error = gr_randtest_not_zero(gr_mat_entry_ptr(rowU, 0, j, ctx), state, ctx);
    if (error != 0) { return error; }
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      if (i >= j) {
        error = gr_set(gr_mat_entry_ptr(L, i, j, ctx), gr_mat_entry_srcptr(colL, i - j, 0, ctx), ctx);
        if (error != 0) { return error; }
      }
      if (j >= i) {
        error = gr_set(gr_mat_entry_ptr(U, i, j, ctx), gr_mat_entry_srcptr(rowU, 0, j - i, ctx), ctx);
        if (error != 0) { return error; }
      }
    }
  }
  error = gr_mat_mul(A, L, U, ctx);
  if (error != 0) { return error; }
  gr_mat_clear(L, ctx);
  gr_mat_clear(U, ctx);
  gr_mat_clear(colL, ctx);
  gr_mat_clear(rowU, ctx);
  return 0;
}

int gr_mat_random_toeplitz(gr_mat_t A, int n, int m, flint_rand_t state, gr_ctx_t ctx) {
  /*Ici je génère de façon non-chalante une matrice toeplitz en prenant 2 vecteur aléatoire que j'utilise pour def la
   * matrice*/
  gr_mat_t col, row;
  int error;
  gr_mat_init(col, n, 1, ctx);
  gr_mat_init(row, 1, m, ctx);

  error = gr_mat_randtest(col, state, ctx);
  if (error != 0) { return error; }
  error = gr_mat_randtest(row, state, ctx);
  if (error != 0) { return error; }
  error = gr_set(gr_mat_entry_ptr(row, 0, 0, ctx), gr_mat_entry_srcptr(col, 0, 0, ctx), ctx);
  if (error != 0) { return error; }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      if (i >= j) {
        error = gr_set(gr_mat_entry_ptr(A, i, j, ctx), gr_mat_entry_srcptr(col, i - j, 0, ctx), ctx);
        if (error != 0) { return error; }
      } else {
        error = gr_set(gr_mat_entry_ptr(A, i, j, ctx), gr_mat_entry_srcptr(row, 0, j - i, ctx), ctx);
        if (error != 0) { return error; }
      }
    }
  }

  gr_mat_clear(col, ctx);
  gr_mat_clear(row, ctx);
  return 0;
}

int gr_mat_quasi_toeplitz_rank(gr_mat_t A, int n, int m, int nb_rand, flint_rand_t state, gr_ctx_t ctx) {
  gr_mat_t temp;
  int error;
  if (nb_rand < 0) { return 1; }

  error = gr_mat_random_quasi_toepitz(A, n, m, state, ctx);
  if (error != 0) { return error; }
  for (int i = 0; i < nb_rand; i++) {
    gr_mat_init(temp, n, m, ctx);
    error = gr_mat_random_quasi_toepitz(temp, n, m, state, ctx);
    if (error != 0) {
      gr_mat_clear(temp, ctx);
      return error;
    }
    error = gr_mat_add(A, A, temp, ctx);
    if (error != 0) {
      gr_mat_clear(temp, ctx);
      return error;
    }
    gr_mat_clear(temp, ctx);
  }

  flint_rand_clear(state);
  return 0;
}
