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
      gr_mul(temp, gr_mat_entry_srcptr(A, i - 1, 0, ctx), gr_mat_entry_srcptr(B, 0, j - 1, ctx), ctx);
      error = gr_sub(gr_mat_entry_ptr(C, i, j, ctx), gr_mat_entry_srcptr(C, i - 1, j - 1, ctx), (gr_srcptr)temp, ctx);
      if (error != 0) { return error; }
    }
  }
  return 0;
}

int gr_multiplication_générateur_déplacement(G_c, H_c, G_a, H_a, G_b, H_b, ctx) {
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

  error = gr_mat_reconstruct_A_safe(G_a, H_a, U, ctx);
  if (error != 0) { return error; }
  error = gr_mat_concat_horizontal(G_c, G_a, U, ctx);
  if (error != 0) { return error; }

  error = gr_mat_reconstruct_A_safe(H_b, G_b, V, ctx);
  if (error != 0) { return error; }
  error = gr_mat_concat_horizontal(H_c, V, H_b, ctx);
  gr_mat_clear(U, ctx);
  gr_mat_clear(V, ctx);
  return error;
}