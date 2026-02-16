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

int gr_multiplication_generateur_deplacement(gr_mat_t G_c, gr_mat_t H_c, gr_mat_t T, gr_mat_t U, gr_mat_t G, gr_mat_t H,
                                             gr_ctx_t ctx) {
  /*
  Donc si j'ai compris cette fois-ci:
  On a 4 générateur de déplacement pour les matrices A et B.
  Pour eviter des multiplication inutile on va juste les utilisé pour avoir en théorie du O(nlog(n))(complexité de FFT
  si je me souviens bien ?) Donc normalement C = G_c * H_c^T et C = A * B donc :
  G_c =  T | W | a , W = A * G ,A = somme( L(gk) * U(Hk))
  H_c = V | H | − b , V = B^t * U ,B = somme( L(hk) * U(gk))
  */
  int error;
  gr_mat_t W, V;

  gr_mat_init(W, gr_mat_nrows(T, ctx), gr_mat_ncols(G, ctx), ctx);
  gr_mat_init(V, gr_mat_nrows(H, ctx), gr_mat_ncols(U, ctx), ctx);
  error = gr_mat_reconstruct_A_safe(W, T, U, ctx);
  if (error != 0) { return error; }
  error = gr_mat_mul(W, W, G, ctx);
  if (error != 0) { return error; }
  error = gr_mat_concat_horizontal(G_c, T, W, ctx);
  if (error != 0) { return error; }

  error = gr_mat_reconstruct_A_safe(V, H, G, ctx);
  if (error != 0) { return error; }
  error = gr_mat_mul(V, V, U, ctx);
  if (error != 0) { return error; }
  error = gr_mat_concat_horizontal(H_c, V, H, ctx);
  gr_mat_clear(W, ctx);
  gr_mat_clear(V, ctx);
  return error;
}

int gr_mat_apply_struct_fast(gr_mat_t Res, gr_mat_t G, gr_mat_t H, gr_mat_t X, gr_ctx_t ctx) {
  slong n = gr_mat_nrows(G, ctx);
  slong alpha = gr_mat_ncols(G, ctx);
  slong m_cols = gr_mat_ncols(X, ctx);

  gr_mat_zero(Res, ctx);

  gr_poly_t poly_v, poly_h, poly_g, tmp_poly, u_v;
  gr_poly_init(poly_v, ctx);
  gr_poly_init(poly_h, ctx);
  gr_poly_init(poly_g, ctx);
  gr_poly_init(tmp_poly, ctx);
  gr_poly_init(u_v, ctx);

  for (slong k = 0; k < alpha; k++) {
    // Redimensionner les polynômes pour accueillir n coefficients
    gr_poly_fit_length(poly_g, n, ctx);
    gr_poly_fit_length(poly_h, n, ctx);
    _gr_poly_set_length(poly_g, n, ctx);
    _gr_poly_set_length(poly_h, n, ctx);

    for (slong i = 0; i < n; i++) {
      // Utilisation de gr_poly_coeff_ptr pour l'accès direct en écriture
      gr_set(gr_poly_coeff_ptr(poly_g, i, ctx), gr_mat_entry_ptr(G, i, k, ctx), ctx);
      gr_set(gr_poly_coeff_ptr(poly_h, i, ctx), gr_mat_entry_ptr(H, n - 1 - i, k, ctx), ctx);
    }

    for (slong j = 0; j < m_cols; j++) {
      gr_poly_fit_length(poly_v, n, ctx);
      _gr_poly_set_length(poly_v, n, ctx);

      for (slong i = 0; i < n; i++) gr_set(gr_poly_coeff_ptr(poly_v, i, ctx), gr_mat_entry_ptr(X, i, j, ctx), ctx);

      // Produit polynomial O(n log n)
      gr_poly_mul(tmp_poly, poly_h, poly_v, ctx);

      // Extraction de la fenêtre pour Toeplitz Supérieure
      gr_poly_fit_length(u_v, n, ctx);
      _gr_poly_set_length(u_v, n, ctx);
      for (slong i = 0; i < n; i++)
        gr_set(gr_poly_coeff_ptr(u_v, i, ctx), gr_poly_coeff_srcptr(tmp_poly, n - 1 + i, ctx), ctx);

      // Produit pour Toeplitz Inférieure
      gr_poly_mul(tmp_poly, poly_g, u_v, ctx);

      // Accumulation
      for (slong i = 0; i < n; i++) {
        gr_add(gr_mat_entry_ptr(Res, i, j, ctx), gr_mat_entry_ptr(Res, i, j, ctx),
               gr_poly_coeff_srcptr(tmp_poly, i, ctx), ctx);
      }
    }
  }

  gr_poly_clear(poly_v, ctx);
  gr_poly_clear(poly_h, ctx);
  gr_poly_clear(poly_g, ctx);
  gr_poly_clear(tmp_poly, ctx);
  gr_poly_clear(u_v, ctx);
  return 0;
}
int gr_mat_apply_Z(gr_mat_t Res, gr_mat_t M, gr_ctx_t ctx) {
  slong n = gr_mat_nrows(M, ctx);
  slong m = gr_mat_ncols(M, ctx);
  gr_mat_zero(Res, ctx);
  for (slong j = 0; j < m; j++) {
    for (slong i = n - 1; i > 0; i--) gr_set(gr_mat_entry_ptr(Res, i, j, ctx), gr_mat_entry_ptr(M, i - 1, j, ctx), ctx);
  }
  return GR_SUCCESS;
}

int gr_mat_apply_Zt(gr_mat_t Res, gr_mat_t M, gr_ctx_t ctx) {
  slong n = gr_mat_nrows(M, ctx);
  slong m = gr_mat_ncols(M, ctx);
  gr_mat_zero(Res, ctx);
  for (slong j = 0; j < m; j++) {
    for (slong i = 0; i < n - 1; i++) gr_set(gr_mat_entry_ptr(Res, i, j, ctx), gr_mat_entry_ptr(M, i + 1, j, ctx), ctx);
  }
  return GR_SUCCESS;
}
int gr_multiplication_generateur_deplacement_fast(gr_mat_t G_c, gr_mat_t H_c, gr_mat_t G_a, gr_mat_t H_a, gr_mat_t G_b,
                                                  gr_mat_t H_b, gr_ctx_t ctx) {
  slong n = gr_mat_nrows(G_a, ctx);
  gr_mat_t W, V, a, b, Tmp, LastCol;

  gr_mat_init(W, n, gr_mat_ncols(G_b, ctx), ctx);
  gr_mat_init(V, n, gr_mat_ncols(H_a, ctx), ctx);
  gr_mat_init(a, n, 1, ctx);
  gr_mat_init(b, n, 1, ctx);
  gr_mat_init(Tmp, n, gr_mat_ncols(G_b, ctx), ctx);
  gr_mat_init(LastCol, n, 1, ctx);

  // 1. W = Z * A * Z^T * G_b
  gr_mat_apply_Zt(Tmp, G_b, ctx);
  gr_mat_apply_struct_fast(W, G_a, H_a, Tmp, ctx);
  gr_mat_apply_Z(W, W, ctx);

  // 2. V = B^T * H_a (Le générateur de B^T est {H_b, G_b})
  gr_mat_apply_struct_fast(V, H_b, G_b, H_a, ctx);

  // 3. Extraction des colonnes de correction (ZA * e_{n-1} et ZB^T * e_{n-1})
  gr_one(gr_mat_entry_ptr(LastCol, n - 1, 0, ctx), ctx);

  gr_mat_apply_struct_fast(a, G_a, H_a, LastCol, ctx);
  gr_mat_apply_Z(a, a, ctx);

  gr_mat_apply_struct_fast(b, H_b, G_b, LastCol, ctx);
  gr_mat_apply_Z(b, b, ctx);

  // 4. Concaténation finale via gr_mat_concat_horizontal
  gr_mat_t G_res, H_res;
  gr_mat_init(G_res, n, gr_mat_ncols(G_a, ctx) + gr_mat_ncols(W, ctx) + 1, ctx);
  gr_mat_init(H_res, n, gr_mat_ncols(V, ctx) + gr_mat_ncols(H_b, ctx) + 1, ctx);

  gr_mat_concat_horizontal(G_res, G_a, W, ctx);
  gr_mat_concat_horizontal(G_res, G_res, a, ctx);

  gr_mat_neg(b, b, ctx);
  gr_mat_concat_horizontal(H_res, V, H_b, ctx);
  gr_mat_concat_horizontal(H_res, H_res, b, ctx);

  gr_mat_set(G_c, G_res, ctx);
  gr_mat_set(H_c, H_res, ctx);

  gr_mat_clear(W, ctx);
  gr_mat_clear(V, ctx);
  gr_mat_clear(a, ctx);
  gr_mat_clear(b, ctx);
  gr_mat_clear(Tmp, ctx);
  gr_mat_clear(LastCol, ctx);
  gr_mat_clear(G_res, ctx);
  gr_mat_clear(H_res, ctx);

  return GR_SUCCESS;
}