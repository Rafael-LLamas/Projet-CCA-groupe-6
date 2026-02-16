#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "matrix_aux.h"
#include "random_toeplitz.h"

#include <stdlib.h>
#include <time.h>

int gr_multiplication_toeplitz(gr_mat_t C, gr_mat_t A, gr_mat_t B, gr_ctx_t ctx) {
  /*
  Si j'ai bien tout compris, on doit d'abord faire sur C le produit scalaire complet pour la ligne 0 et la colonne 0.
  Puis je prends en compte que la 1ere colonne et ligne  de A et B.
  Avec ca je dois suivre la formule : C i,j = Ci-1,j-1 + (Ai,n-1 x Bn-1,j) - (Ai-1,0 x B0,j-1), sachant que les matrices
  sont Toeplitz donc je peux simplifier A et B par leurs vecteur et juste prendre la bonne position
  */
  int error;
  for (int i = 0; i < gr_mat_nrows(C, ctx); i++) {
    for (int j = 0; j < gr_mat_ncols(A, ctx); j++) {
      error = gr_mul(gr_mat_entry_ptr(C, i, 0, ctx), gr_mat_entry_srcptr(A, i, j, ctx),
                     gr_mat_entry_srcptr(B, j, 0, ctx), ctx);
      if (error != 0) { return error; }
    }
  }
  for (int j = 1; j < gr_mat_ncols(C, ctx); j++) {
    for (int i = 0; i < gr_mat_ncols(A, ctx); i++) {
      error = gr_mul(gr_mat_entry_ptr(C, 0, j, ctx), gr_mat_entry_srcptr(A, 0, i, ctx),
                     gr_mat_entry_srcptr(B, i, j, ctx), ctx);
      if (error != 0) { return error; }
    }
  }
  gr_ptr temp;
  GR_TMP_INIT(temp, ctx);
  for (int i = 1; i < gr_mat_nrows(C, ctx); i++) {
    for (int j = 1; j < gr_mat_ncols(C, ctx); j++) {
      error = gr_mul(gr_mat_entry_ptr(C, i, j, ctx), gr_mat_entry_srcptr(A, i, gr_mat_ncols(A, ctx) - 1, ctx),
                     gr_mat_entry_srcptr(B, gr_mat_ncols(A, ctx) - 1, j, ctx), ctx);
      if (error != 0) { return error; }
      error = gr_add(gr_mat_entry_ptr(C, i, j, ctx), gr_mat_entry_srcptr(C, i - 1, j - 1, ctx),
                     gr_mat_entry_ptr(C, i, j, ctx), ctx);
      if (error != 0) { return error; }
      FLINT_CHECK(gr_mul(temp, gr_mat_entry_srcptr(A, i - 1, 0, ctx), gr_mat_entry_srcptr(B, 0, j - 1, ctx), ctx));
      error = gr_sub(gr_mat_entry_ptr(C, i, j, ctx), gr_mat_entry_srcptr(C, i - 1, j - 1, ctx), (gr_srcptr)temp, ctx);
      if (error != 0) { return error; }
    }
  }
  return 0;
}

int gr_multiplication_generateur_deplacement(gr_mat_t G_c, gr_mat_t H_c, gr_mat_t G_a, gr_mat_t H_a, gr_mat_t G_b,
                                             gr_mat_t H_b, gr_ctx_t ctx) {
  /*
  Donc si j'ai compris cette fois-ci:
  On a 4 générateur de déplacement pour les matrices A et B.
  Pour eviter des multiplication inutile on va juste les utilisé pour avoir en théorie du O(nlog(n))(complexité de FFT
  si je me souviens bien ?) Donc normalement C = G_c * H_c^T et C = A * B donc :
  G_c = G_A | U , U = A * G_b ,A = somme( L(gk) * U(Hk))
  H_c = V | H_b , V = B^t * H_a ,B = somme( L(hk) * U(gk))
  */
  int error;
  gr_mat_t U, V;

  gr_mat_init(U, gr_mat_nrows(G_a, ctx), gr_mat_ncols(G_b, ctx), ctx);
  gr_mat_init(V, gr_mat_nrows(H_b, ctx), gr_mat_ncols(H_a, ctx), ctx);
  error = gr_mat_reconstruct_A_safe(U, G_a, H_a, ctx);
  if (error != 0) { return error; }
  error = gr_mat_mul(U, U, G_b, ctx);
  if (error != 0) { return error; }
  error = gr_mat_concat_horizontal(G_c, G_a, U, ctx);
  if (error != 0) { return error; }

  error = gr_mat_reconstruct_A_safe(V, H_b, G_b, ctx);
  if (error != 0) { return error; }
  error = gr_mat_mul(V, V, H_a, ctx);
  if (error != 0) { return error; }
  error = gr_mat_concat_horizontal(H_c, V, H_b, ctx);
  gr_mat_clear(U, ctx);
  gr_mat_clear(V, ctx);
  return error;
}

int test_multiplication_generateurs() {
  int error;
  gr_mat_t A, B, C, G_a, H_a, G_b, H_b, G_c, H_c;
  gr_ctx_t ctx;
  flint_rand_t state;
  gr_ctx_init_nmod(ctx, 47);
  flint_rand_init(state);
  error = random_toeplitz(A, gr_mat_nrows(A, ctx), gr_mat_ncols(A, ctx), state, ctx);
  if (error != 0) { return error; }
  error = random_toeplitz(B, gr_mat_nrows(B, ctx), gr_mat_ncols(B, ctx), state, ctx);
  if (error != 0) { return error; }
  gr_mat_init(C, gr_mat_nrows(A, ctx), gr_mat_ncols(B, ctx), ctx);
  error = gr_mat_G_H(G_a, H_a, A, ctx);
  if (error != 0) { return error; }
  error = gr_mat_G_H(G_b, H_b, B, ctx);
  if (error != 0) { return error; }
  error = gr_multiplication_toeplitz(C, A, B, ctx);
  if (error != 0) { return error; }
  flint_printf(":------------------------Matrices C faite sans les générateurs---------------------------------:");
  flint_printf("Matrice C = \n");
  gr_mat_print(C, ctx);
  flint_printf("\n");
  gr_mat_G_H(C, G_c, H_c, ctx);
  flint_printf("Matrice G_c = \n");
  gr_mat_print(G_c, ctx);
  flint_printf("\n");
  flint_printf("Matrice H_c = \n");
  gr_mat_print(H_c, ctx);
  flint_printf("\n");
  gr_multiplication_generateur_deplacement(G_c, H_c, G_a, H_a, G_b, H_b, ctx);
  flint_printf(":------------------------Matrices C faite avec les générateurs---------------------------------:");
  flint_printf("Matrice G_c = \n");
  gr_mat_print(G_c, ctx);
  flint_printf("\n");
  flint_printf("Matrice H_c = \n");
  gr_mat_print(H_c, ctx);
  flint_printf("\n");
  gr_mat_reconstruct_A_safe(C, G_c, H_c, ctx);
  flint_printf("Matrice C = \n");
  gr_mat_print(C, ctx);
  flint_printf("\n");

  return error;
}