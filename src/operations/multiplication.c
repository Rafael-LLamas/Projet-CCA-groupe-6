#include "matrix_aux.h"

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include <stdlib.h>
#include <time.h>

int gr_mat_mul_vector(gr_mat_t Res, gr_mat_t G, gr_mat_t H, gr_mat_t X, gr_ctx_t ctx) {
  int error;
  slong rows_G = gr_mat_nrows(G, ctx);
  slong rows_X = gr_mat_nrows(X, ctx);
  slong alpha = gr_mat_ncols(G, ctx);
  slong m_cols = gr_mat_ncols(X, ctx);

  // Res doit être de taille (rows_G x m_cols)
  error = gr_mat_zero(Res, ctx);
  if (error != 0) return error;

  gr_poly_t poly_v, poly_h, poly_g, tmp_poly, u_v;
  gr_poly_init(poly_v, ctx);
  gr_poly_init(poly_h, ctx);
  gr_poly_init(poly_g, ctx);
  gr_poly_init(tmp_poly, ctx);
  gr_poly_init(u_v, ctx);

  for (slong k = 0; k < alpha; k++) {
    // Préparation des polynômes du générateur (taille rows_G)
    gr_poly_fit_length(poly_g, rows_G, ctx);
    gr_poly_fit_length(poly_h, rows_G, ctx);
    _gr_poly_set_length(poly_g, rows_G, ctx);
    _gr_poly_set_length(poly_h, rows_G, ctx);

    for (slong i = 0; i < rows_G; i++) {
      error = gr_set(gr_poly_coeff_ptr(poly_g, i, ctx), gr_mat_entry_srcptr(G, i, k, ctx), ctx);
      if (error != 0) return error;
      // On utilise rows_G pour l'inversion d'indice de H
      error = gr_set(gr_poly_coeff_ptr(poly_h, i, ctx), gr_mat_entry_srcptr(H, rows_G - 1 - i, k, ctx), ctx);
      if (error != 0) return error;
    }

    for (slong j = 0; j < m_cols; j++) {
      // Préparation du polynôme vecteur (taille rows_X)
      gr_poly_fit_length(poly_v, rows_X, ctx);
      _gr_poly_set_length(poly_v, rows_X, ctx);

      for (slong i = 0; i < rows_X; i++) {
        error = gr_set(gr_poly_coeff_ptr(poly_v, i, ctx), gr_mat_entry_srcptr(X, i, j, ctx), ctx);
        if (error != 0) return error;
      }

      // Convolution pour la partie Toeplitz
      error = gr_poly_mul(tmp_poly, poly_h, poly_v, ctx);
      if (error != 0) return error;

      // Extraction de la fenêtre (doit être cohérente avec rows_X)
      gr_poly_fit_length(u_v, rows_X, ctx);
      _gr_poly_set_length(u_v, rows_X, ctx);
      for (slong i = 0; i < rows_X; i++) {
        // n-1 devient rows_X-1 pour gérer les matrices non-carrées
        error = gr_set(gr_poly_coeff_ptr(u_v, i, ctx), gr_poly_coeff_srcptr(tmp_poly, rows_X - 1 + i, ctx), ctx);
        if (error != 0) return error;
      }

      error = gr_poly_mul(tmp_poly, poly_g, u_v, ctx);
      if (error != 0) return error;

      // Accumulation dans le résultat (taille rows_G)
      for (slong i = 0; i < rows_G; i++) {
        if (i < gr_poly_length(tmp_poly, ctx)) {
          error = gr_add(gr_mat_entry_ptr(Res, i, j, ctx), gr_mat_entry_srcptr(Res, i, j, ctx),
                         gr_poly_coeff_srcptr(tmp_poly, i, ctx), ctx);
          if (error != 0) return error;
        }
      }
    }
  }

  gr_poly_clear(poly_v, ctx);
  gr_poly_clear(poly_h, ctx);
  gr_poly_clear(poly_g, ctx);
  gr_poly_clear(tmp_poly, ctx);
  gr_poly_clear(u_v, ctx);
  return GR_SUCCESS;
}

int gr_mat_mul_generator(gr_mat_t G_c, gr_mat_t H_c, gr_mat_t G_a, gr_mat_t H_a, gr_mat_t G_b, gr_mat_t H_b,
                         gr_ctx_t ctx) {
  slong n = gr_mat_nrows(G_a, ctx); // Lignes de A
  slong m = gr_mat_nrows(G_b, ctx); // Colonnes de A
  slong k = gr_mat_nrows(H_b, ctx); // colonnes de B

  gr_mat_t W, V, a, b, Tmp, LastCol_m, LastCol_k;
  int error;

  // Initialisation avec les dimensions respectives
  gr_mat_init(W, n, gr_mat_ncols(G_b, ctx), ctx);
  gr_mat_init(V, k, gr_mat_ncols(H_a, ctx), ctx);
  gr_mat_init(a, n, 1, ctx);
  gr_mat_init(b, k, 1, ctx);
  gr_mat_init(Tmp, m, gr_mat_ncols(G_b, ctx), ctx);
  gr_mat_init(LastCol_m, m, 1, ctx);
  gr_mat_init(LastCol_k, k, 1, ctx);

  // 1. W = Z * A * Z^T * G_b
  error = gr_mat_apply_Zt(Tmp, G_b, ctx);
  if (error != 0) { return error; }
  error = gr_mat_mul_vector(W, G_a, H_a, Tmp, ctx);
  if (error != 0) { return error; }
  error = gr_mat_apply_Z(W, W, ctx);
  if (error != 0) { return error; }

  // 2. V = B^T * H_a (Le générateur de B^T est {H_b, G_b})
  error = gr_mat_mul_vector(V, H_b, G_b, H_a, ctx);
  if (error != 0) { return error; }

  // 3. Correction colonnes
  // Vecteur e_{m-1} pour la dimension de A
  error = gr_one(gr_mat_entry_ptr(LastCol_m, m - 1, 0, ctx), ctx);
  if (error != 0) { return error; }
  error = gr_mat_mul_vector(a, G_a, H_a, LastCol_m, ctx);
  if (error != 0) { return error; }
  error = gr_mat_apply_Z(a, a, ctx);
  if (error != 0) { return error; }

  // Vecteur e_{k-1} pour la dimension de B^T
  error = gr_one(gr_mat_entry_ptr(LastCol_k, k - 1, 0, ctx), ctx);
  if (error != 0) { return error; }
  error = gr_mat_mul_vector(b, H_b, G_b, LastCol_k, ctx);
  if (error != 0) { return error; }
  error = gr_mat_apply_Z(b, b, ctx);
  if (error != 0) { return error; }

  // 4. Concaténation (G_c: n x (ra+rb+1), H_c: k x (ra+rb+1))
  gr_mat_t G_res, H_res, G_temp, H_temp;
  gr_mat_init(G_res, n, gr_mat_ncols(G_a, ctx) + gr_mat_ncols(W, ctx) + 1, ctx);
  gr_mat_init(G_temp, n, gr_mat_ncols(G_a, ctx) + gr_mat_ncols(W, ctx), ctx);
  gr_mat_init(H_res, k, gr_mat_ncols(V, ctx) + gr_mat_ncols(H_b, ctx) + 1, ctx);
  gr_mat_init(H_temp, k, gr_mat_ncols(V, ctx) + gr_mat_ncols(H_b, ctx), ctx);

  error = gr_mat_concat_horizontal(G_temp, G_a, W, ctx);
  if (error != 0) { return error; }
  error = gr_mat_concat_horizontal(G_res, G_temp, a, ctx);
  if (error != 0) { return error; }

  error = gr_mat_neg(b, b, ctx);
  if (error != 0) { return error; }
  error = gr_mat_concat_horizontal(H_temp, V, H_b, ctx);
  if (error != 0) { return error; }
  error = gr_mat_concat_horizontal(H_res, H_temp, b, ctx);
  if (error != 0) { return error; }

  // On remplace les anciennes matrices par les nouvelles
  error = gr_mat_set(G_c, G_res, ctx);
  if (error != 0) return error;
  error = gr_mat_set(H_c, H_res, ctx);
  if (error != 0) return error;

  gr_mat_clear(W, ctx);
  gr_mat_clear(V, ctx);
  gr_mat_clear(a, ctx);
  gr_mat_clear(b, ctx);
  gr_mat_clear(Tmp, ctx);
  gr_mat_clear(LastCol_m, ctx);
  gr_mat_clear(LastCol_k, ctx);
  gr_mat_clear(G_temp, ctx);
  gr_mat_clear(H_temp, ctx);
  gr_mat_clear(G_res, ctx);
  gr_mat_clear(H_res, ctx);

  return GR_SUCCESS;
}