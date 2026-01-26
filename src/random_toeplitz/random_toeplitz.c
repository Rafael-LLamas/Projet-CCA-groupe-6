#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include <stdlib.h>

void random_toeplitz(gr_mat_t A, int n, int m, gr_ctx_t ctx) {
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
  // Génération d'une matrice random et de Z a zero
  flint_rand_t state;
  gr_mat_t Z, ZT, LU;
  slong *rank, *P;
  gr_mat_init(Z, n, m, ctx);
  gr_mat_init(LU, n, m, ctx);
  int error = gr_mat_zero(Z, ctx);
  flint_rand_init(state);
  gr_mat_init(A, 10, 10, ctx);
  error = gr_mat_randtest(A, state, ctx);

  if (error) { // erreur ??
    gr_mat_clear(Z, ctx);
    gr_mat_clear(A, ctx);
    gr_ctx_clear(ctx);
    flint_rand_clear(state);
    exit(error);
  }
  // Décomposition LU
  error = gr_mat_lu(rank, P, LU, A, 0, ctx);
  if (error) { // erreur ??
    gr_mat_clear(Z, ctx);
    gr_mat_clear(A, ctx);
    gr_ctx_clear(ctx);
    flint_rand_clear(state);
    exit(error);
  }
  // clear l'inutiliser
  flint_rand_clear(state);
  gr_mat_clear(Z, ctx);
}

int main() {
  gr_ctx_t ctx;
  gr_ctx_init_fmpz(ctx);
  gr_mat_t ran;
  random_toeplitz(ran, 10, 10, ctx);
  flint_printf("mat1 = ");
  gr_mat_print(ran, ctx);
  flint_printf("\n");
  gr_mat_clear(ran, ctx);
  gr_ctx_clear(ctx);
}