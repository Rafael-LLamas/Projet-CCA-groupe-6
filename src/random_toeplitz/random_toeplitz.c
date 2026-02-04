#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include <stdlib.h>

int random_quasi_toeplitz(gr_mat_t A, int n, int m, flint_rand_t state, gr_ctx_t ctx) {
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

int random_toeplitz(gr_mat_t A, int n, int m, flint_rand_t state, gr_ctx_t ctx) {
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
  error = gr_set(gr_mat_entry_ptr(row, 0, 0, ctx), gr_mat_entry_ptr(col, 0, 0, ctx), ctx);
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
int rand_quasi_toeplitz(gr_mat_t A, int n, int m, int nb_rand, gr_ctx_t ctx) {
  /*
  Pour palié le manque de random je vais générer plusieurs matrice quasi-toeplitz et les additionner entre elles,
  nh'ésitez pas a faire varier nb_rand.
  !!!Attention, je pense qu'il y a une possiblilité que la matrice ne soit plus quasi-toeplitz si nb_rand est trop
  grand!!!(faire une test pour trouver la limite ??)
  */
  flint_rand_t state;
  gr_mat_t temp;
  int error;
  if (nb_rand < 0) { return 1; }
  flint_rand_init(state);

  error = random_quasi_toeplitz(A, n, m, state, ctx);
  if (error != 0) { return error; }
  for (int i = 0; i < nb_rand; i++) {
    gr_mat_init(temp, n, m, ctx);
    error = random_quasi_toeplitz(temp, n, m, state, ctx);
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

int test_random_toeplitz() {
  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, 47);
  int error;
  gr_mat_t ran, D;
  slong *rank = flint_malloc(sizeof(slong));
  gr_mat_init(ran, 5, 5, ctx);
  error = rand_quasi_toeplitz(ran, 5, 5, 1, ctx);
  if (error != 0) {
    gr_mat_clear(ran, ctx);
    gr_ctx_clear(ctx);
    flint_free(rank);
    exit(error);
  }
  error = gr_mat_rank(rank, ran, ctx);
  if (error != 0) {
    flint_printf("%d\n", error);
    gr_mat_clear(ran, ctx);
    gr_ctx_clear(ctx);
    flint_free(rank);
    exit(error);
  }
  flint_printf("----------------------------------------------\nResultat quasi = ");
  gr_mat_print(ran, ctx);
  flint_printf("\n");
  flint_printf("rank = %wd \n", *rank);
  gr_mat_clear(ran, ctx);
  gr_mat_init(ran, 5, 5, ctx);
  error = rand_quasi_toeplitz(ran, 5, 5, 4, ctx);
  if (error != 0) {
    gr_mat_clear(ran, ctx);
    gr_ctx_clear(ctx);
    flint_free(rank);
    exit(error);
  }

  flint_printf("----------------------------------------------\nResultat quasi = ");
  gr_mat_print(ran, ctx);
  flint_printf("\n");
  gr_mat_clear(ran, ctx);
  gr_mat_init(ran, 5, 5, ctx);
  error = rand_quasi_toeplitz(ran, 5, 5, 2, ctx);
  if (error != 0) {
    gr_mat_clear(ran, ctx);
    gr_ctx_clear(ctx);
    flint_free(rank);
    exit(error);
  }
  flint_printf("----------------------------------------------\nResultat quasi = ");
  gr_mat_print(ran, ctx);
  flint_printf("\n");
  gr_mat_clear(ran, ctx);
  gr_mat_init(ran, 5, 5, ctx);
  error = rand_quasi_toeplitz(ran, 5, 5, 3, ctx);
  if (error != 0) {
    gr_mat_clear(ran, ctx);
    gr_ctx_clear(ctx);
    flint_free(rank);
    exit(error);
  }

  flint_printf("----------------------------------------------\nResultat quasi = ");
  gr_mat_print(ran, ctx);
  flint_printf("\n");
  gr_mat_clear(ran, ctx);
  gr_mat_init(ran, 15, 15, ctx);
  error = rand_quasi_toeplitz(ran, 15, 15, 0, ctx);
  if (error != 0) {
    gr_mat_clear(ran, ctx);
    gr_ctx_clear(ctx);
    flint_free(rank);
    exit(error);
  }
  gr_mat_init(D, 15, 15, ctx);
  error = gr_mat_displacement_square_safe(D, ran, ctx);
  if (error != 0) {
    gr_mat_clear(ran, ctx);
    gr_mat_clear(D, ctx);
    gr_ctx_clear(ctx);
    flint_free(rank);
    exit(error);
  }
  error = gr_mat_rank(rank, D, ctx);
  if (error != 0) {
    gr_mat_clear(ran, ctx);
    gr_mat_clear(D, ctx);
    gr_ctx_clear(ctx);
    flint_free(rank);
    exit(error);
  }
  flint_printf("----------------------------------------------\nResultat quasi = ");
  gr_mat_print(ran, ctx);
  flint_printf("\n\n\n");
  flint_printf(" deplacement =  ");
  gr_mat_print(D, ctx);
  flint_printf("\n");
  flint_printf("rank de déplcement = %wd \n", *rank);
  gr_mat_clear(D, ctx);
  gr_mat_clear(ran, ctx);
  gr_ctx_clear(ctx);
  flint_free(rank);
  return GR_SUCCESS;
}