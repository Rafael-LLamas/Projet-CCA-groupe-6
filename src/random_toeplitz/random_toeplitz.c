#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include <stdlib.h>

int gr_mat_random_quasi_toepitz(gr_mat_t A, flint_rand_t state, gr_ctx_t ctx) {
  /*
  Mon sommeil va en patir ......
  Donc apres avoir compris mon problème ici je génère des quasi toeplitz en suivant la logique que mets matrices auront
  un rang faible. Pour ca je génère 2 matrice toeplitz, L et U puis je l'ai multiplie entre elle pour avoir A.
  */

  gr_mat_t L, U, colL, rowU;
  int error;
  gr_mat_init(L, gr_mat_nrows(A, ctx), gr_mat_ncols(A, ctx), ctx);
  gr_mat_init(U, gr_mat_ncols(A, ctx), gr_mat_ncols(A, ctx), ctx);
  gr_mat_init(colL, gr_mat_nrows(A, ctx), 1, ctx);
  gr_mat_init(rowU, 1, gr_mat_ncols(A, ctx), ctx);
  error = gr_mat_zero(L, ctx);
  if (error != GR_SUCCESS) {
    gr_mat_clear(L, ctx);
    gr_mat_clear(U, ctx);
    gr_mat_clear(colL, ctx);
    gr_mat_clear(rowU, ctx);
    return error;
  }
  error = gr_mat_zero(U, ctx);
  if (error != GR_SUCCESS) {
    gr_mat_clear(L, ctx);
    gr_mat_clear(U, ctx);
    gr_mat_clear(colL, ctx);
    gr_mat_clear(rowU, ctx);
    return error;
  }

  for (int i = 0; i < gr_mat_nrows(A, ctx); i++) {
    error = gr_randtest_not_zero(gr_mat_entry_ptr(colL, i, 0, ctx), state, ctx);
    if (error != GR_SUCCESS) {
      gr_mat_clear(L, ctx);
      gr_mat_clear(U, ctx);
      gr_mat_clear(colL, ctx);
      gr_mat_clear(rowU, ctx);
      return error;
    }
  }
  for (int j = 0; j < gr_mat_ncols(A, ctx); j++) {
    error = gr_randtest_not_zero(gr_mat_entry_ptr(rowU, 0, j, ctx), state, ctx);
    if (error != GR_SUCCESS) {
      gr_mat_clear(L, ctx);
      gr_mat_clear(U, ctx);
      gr_mat_clear(colL, ctx);
      gr_mat_clear(rowU, ctx);
      return error;
    }
  }

  for (int i = 0; i < gr_mat_nrows(A, ctx); i++) {
    for (int j = 0; j < gr_mat_ncols(A, ctx); j++) {
      if (i >= j) {
        error = gr_set(gr_mat_entry_ptr(L, i, j, ctx), gr_mat_entry_srcptr(colL, i - j, 0, ctx), ctx);
        if (error != GR_SUCCESS) {
          gr_mat_clear(L, ctx);
          gr_mat_clear(U, ctx);
          gr_mat_clear(colL, ctx);
          gr_mat_clear(rowU, ctx);
          return error;
        }
      }
      if (j >= i) {
        error = gr_set(gr_mat_entry_ptr(U, i, j, ctx), gr_mat_entry_srcptr(rowU, 0, j - i, ctx), ctx);
        if (error != GR_SUCCESS) {
          gr_mat_clear(L, ctx);
          gr_mat_clear(U, ctx);
          gr_mat_clear(colL, ctx);
          gr_mat_clear(rowU, ctx);
          return error;
        }
      }
    }
  }
  error = gr_mat_mul(A, L, U, ctx);
  if (error != GR_SUCCESS) { return error; }
  gr_mat_clear(L, ctx);
  gr_mat_clear(U, ctx);
  gr_mat_clear(colL, ctx);
  gr_mat_clear(rowU, ctx);
  return 0;
}

int gr_mat_random_toeplitz(gr_mat_t A, flint_rand_t state, gr_ctx_t ctx) {
  /*Ici je génère de façon non-chalante une matrice toeplitz en prenant 2 vecteur aléatoire que j'utilise pour def la
   * matrice*/
  gr_mat_t col, row;
  int error;
  gr_mat_init(col, gr_mat_nrows(A, ctx), 1, ctx);
  gr_mat_init(row, 1, gr_mat_ncols(A, ctx), ctx);

  for (int i = 0; i < gr_mat_nrows(A, ctx); i++) {
    error = gr_randtest_not_zero(gr_mat_entry_ptr(col, i, 0, ctx), state, ctx);
    if (error != GR_SUCCESS) {
      gr_mat_clear(col, ctx);
      gr_mat_clear(row, ctx);
      return error;
    }
  }
  for (int j = 1; j < gr_mat_ncols(A, ctx); j++) {
    error = gr_randtest_not_zero(gr_mat_entry_ptr(row, 0, j, ctx), state, ctx);
    if (error != GR_SUCCESS) {
      gr_mat_clear(col, ctx);
      gr_mat_clear(row, ctx);
      return error;
    }
  }
  if (error != GR_SUCCESS) {
    gr_mat_clear(col, ctx);
    gr_mat_clear(row, ctx);
    return error;
  }
  for (int i = 0; i < gr_mat_nrows(A, ctx); i++) {
    for (int j = 0; j < gr_mat_ncols(A, ctx); j++) {
      if (i >= j) {
        error = gr_set(gr_mat_entry_ptr(A, i, j, ctx), gr_mat_entry_srcptr(col, i - j, 0, ctx), ctx);
        if (error != GR_SUCCESS) {
          gr_mat_clear(col, ctx);
          gr_mat_clear(row, ctx);
          return error;
        }
      } else {
        error = gr_set(gr_mat_entry_ptr(A, i, j, ctx), gr_mat_entry_srcptr(row, 0, j - i, ctx), ctx);
        if (error != GR_SUCCESS) {
          gr_mat_clear(col, ctx);
          gr_mat_clear(row, ctx);
          return error;
        }
      }
    }
  }

  gr_mat_clear(col, ctx);
  gr_mat_clear(row, ctx);
  return 0;
}

int gr_mat_quasi_toeplitz_rank(gr_mat_t A, int nb_rand, flint_rand_t state, gr_ctx_t ctx) {
  gr_mat_t temp;
  int error;
  if (nb_rand < 0) { return 1; }

  error = gr_mat_random_quasi_toepitz(A, state, ctx);
  if (error != GR_SUCCESS) { return error; }
  for (int i = 1; i < nb_rand; i++) {
    gr_mat_init(temp, gr_mat_nrows(A, ctx), gr_mat_ncols(A, ctx), ctx);
    error = gr_mat_random_quasi_toepitz(temp, state, ctx);
    if (error != GR_SUCCESS) {
      gr_mat_clear(temp, ctx);
      return error;
    }
    error = gr_mat_add(A, A, temp, ctx);
    if (error != GR_SUCCESS) {
      gr_mat_clear(temp, ctx);
      return error;
    }
    gr_mat_clear(temp, ctx);
  }

  return 0;
}

int gr_mat_random_generator_toeplitz(gr_mat_t G, gr_mat_t H, flint_rand_t state, gr_ctx_t ctx) {
  int error = GR_SUCCESS;
  for (int i = 0; i < gr_mat_nrows(G, ctx); i++) {
    error = gr_randtest_not_zero(gr_mat_entry_ptr(G, i, 0, ctx), state, ctx);
    if (error != GR_SUCCESS) { return error; }
  }
  for (int j = 1; j < gr_mat_nrows(H, ctx); j++) {
    error = gr_randtest_not_zero(gr_mat_entry_ptr(H, j, 0, ctx), state, ctx);
    if (error != GR_SUCCESS) { return error; }
  }
  error = gr_set(gr_mat_entry_ptr(H, 0, 0, ctx), gr_mat_entry_srcptr(G, 0, 0, ctx), ctx);
  return error;
}
