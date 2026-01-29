#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include <stdlib.h>

void random_toeplitz(gr_mat_t A, int n, int m, flint_rand_t state, gr_ctx_t ctx) {
  /*
  Pour generer une matrice toeplitz il me faut :
  1- generer une matrice aléatoire de taille n*m
  2- faire une décomposition LU
  3-Vérifier la décomposition
  4-Créer Z
  5- faire A - ZAZ^t
  6-Vérifier le rank
  Normalement je vais obtenir une matrice quasi-toeplitz (normalement)
  */
  // Génération d'une matrice random

  gr_mat_t Z, LU;
  slong *rank, *P;
  gr_mat_init(A, n, m, ctx);
  P = flint_malloc(n * sizeof(slong));
  rank = flint_malloc(sizeof(slong));
  int error = gr_mat_randtest(A, state, ctx);
  flint_printf("----------------------------------------------\nA = ");
  gr_mat_print(A, ctx);
  flint_printf("\n");
  if (error != 0) { // erreur ??
    gr_mat_clear(Z, ctx);
    gr_mat_clear(A, ctx);
    gr_ctx_clear(ctx);
    flint_rand_clear(state);
    flint_free(P);
    flint_free(rank);
    exit(error);
  }
  // init LU et Z
  error = gr_mat_init_set(LU, A, ctx);
  gr_mat_init(Z, n, m, ctx);
  error = gr_mat_zero(Z, ctx);

  int i = 1, j = 0;

  while (i < n && j < m - 1) {
    error = gr_set_si(gr_mat_entry_ptr(Z, i, j, ctx), 1, ctx);
    i++;
    j++;
  }

  flint_printf("----------------------------------------------\nZ = ");
  gr_mat_print(Z, ctx);
  flint_printf("\n");

  // Décomposition LU
  error = gr_mat_displacement(A, A, ctx);
  error = gr_mat_lu(rank, P, A, A, 0, ctx);
  flint_printf("----------------------------------------------\nLU = ");
  gr_mat_print(LU, ctx);
  flint_printf("\nrank LU = %{slong}", *rank);
  flint_printf("\n");

  if (error != 0) { // erreur ??
    gr_mat_clear(Z, ctx);
    flint_printf("%d\n", error);
    gr_mat_clear(A, ctx);
    gr_ctx_clear(ctx);
    flint_rand_clear(state);
    flint_free(P);
    flint_free(rank);
    exit(error);
  }

  // // création de la toeplitz like avec LU-ZLUZ^T
  // error = gr_mat_mul(A, LU, Z, ctx); // A = ZLU
  // if (error != 0) {                  // erreur ??
  //   gr_mat_clear(Z, ctx);
  //   flint_printf("%d\n", error);
  //   gr_mat_clear(A, ctx);
  //   gr_ctx_clear(ctx);
  //   flint_rand_clear(state);
  //   flint_free(P);
  //   flint_free(rank);
  //   exit(error);
  // }
  // flint_printf("----------------------------------------------\nZLU = ");
  // gr_mat_print(A, ctx);
  // flint_printf("\n");
  // error = gr_mat_transpose(Z, Z, ctx); // Z == Z^T
  // if (error != 0) {                    // erreur ??
  //   gr_mat_clear(Z, ctx);
  //   flint_printf("%d\n", error);
  //   gr_mat_clear(A, ctx);
  //   gr_ctx_clear(ctx);
  //   flint_rand_clear(state);
  //   flint_free(P);
  //   flint_free(rank);
  //   exit(error);
  // }
  // flint_printf("----------------------------------------------\nZ^T = ");
  // gr_mat_print(Z, ctx);
  // flint_printf("\n");

  // error = gr_mat_mul(A, A, Z, ctx); // A = ZLUZ^T
  // if (error != 0) {                 // erreur ??
  //   gr_mat_clear(Z, ctx);
  //   flint_printf("%d\n", error);
  //   gr_mat_clear(A, ctx);
  //   gr_ctx_clear(ctx);
  //   flint_rand_clear(state);
  //   flint_free(P);
  //   flint_free(rank);
  //   exit(error);
  // }
  // flint_printf("----------------------------------------------\nZLUZ^T = ");
  // gr_mat_print(A, ctx);
  // error = gr_mat_rank(rank, A, ctx);
  // flint_printf("\nrank ZLUZ^T = %{slong}", *rank);
  // flint_printf("\n");
  // error = gr_mat_sub(A, LU, A, ctx); // A = LU-ZLUZ^T
  // if (error != 0) {                  // erreur ??

  //   gr_mat_clear(Z, ctx);
  //   flint_printf("%d\n", error);
  //   gr_mat_clear(A, ctx);
  //   gr_ctx_clear(ctx);
  //   flint_rand_clear(state);
  //   flint_free(P);
  //   flint_free(rank);
  //   exit(error);
  // }

  // clear l'inutiliser

  gr_mat_clear(Z, ctx);
  gr_mat_clear(LU, ctx);
  flint_free(P);
  flint_free(rank);
  flint_rand_clear(state);
}

int test_random_toeplitz() {
  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, 1009);
  gr_mat_t ran;
  flint_rand_t state;
  flint_rand_init(state);
  random_toeplitz(ran, 5, 5, state, ctx);
  flint_printf("----------------------------------------------\nResultat = ");
  gr_mat_print(ran, ctx);
  flint_printf("\n");
  gr_mat_clear(ran, ctx);
  gr_ctx_clear(ctx);
  return GR_SUCCESS;
}