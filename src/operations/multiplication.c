#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"

#include <stdlib.h>
#include <time.h>

int gr_mat_mul_vector(gr_mat_t Res, gr_mat_t G, gr_mat_t H, gr_mat_t X, gr_ctx_t ctx) {
  int error;
  slong n = gr_mat_nrows(G, ctx);
  slong alpha = gr_mat_ncols(G, ctx);
  slong m_cols = gr_mat_ncols(X, ctx);

  error = gr_mat_zero(Res, ctx);
  if (error != 0) { return error; }

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
      error = gr_set(gr_poly_coeff_ptr(poly_g, i, ctx), gr_mat_entry_srcptr(G, i, k, ctx), ctx);
      if (error != 0) { return error; }
      error = gr_set(gr_poly_coeff_ptr(poly_h, i, ctx), gr_mat_entry_srcptr(H, n - 1 - i, k, ctx), ctx);
      if (error != 0) { return error; }
    }

    for (slong j = 0; j < m_cols; j++) {
      gr_poly_fit_length(poly_v, n, ctx);
      _gr_poly_set_length(poly_v, n, ctx);

      for (slong i = 0; i < n; i++) {
        error = gr_set(gr_poly_coeff_ptr(poly_v, i, ctx), gr_mat_entry_srcptr(X, i, j, ctx), ctx);
        if (error != 0) { return error; }
      }

      // Produit polynomial O(n log n)
      error = gr_poly_mul(tmp_poly, poly_h, poly_v, ctx);
      if (error != 0) { return error; }

      // Extraction de la fenêtre pour Toeplitz Supérieure
      gr_poly_fit_length(u_v, n, ctx);
      _gr_poly_set_length(u_v, n, ctx);
      for (slong i = 0; i < n; i++) {
        error = gr_set(gr_poly_coeff_ptr(u_v, i, ctx), gr_poly_coeff_srcptr(tmp_poly, n - 1 + i, ctx), ctx);
        if (error != 0) { return error; }
      }

      // Produit pour Toeplitz Inférieure
      error = gr_poly_mul(tmp_poly, poly_g, u_v, ctx);
      if (error != 0) { return error; }

      // Accumulation
      for (slong i = 0; i < n; i++) {
        error = gr_add(gr_mat_entry_ptr(Res, i, j, ctx), gr_mat_entry_srcptr(Res, i, j, ctx),
                       gr_poly_coeff_srcptr(tmp_poly, i, ctx), ctx);
        if (error != 0) { return error; }
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
  int error;
  slong n = gr_mat_nrows(M, ctx);
  slong m = gr_mat_ncols(M, ctx);

  for (slong j = 0; j < m; j++) {
    // On part du bas vers le haut pour ne pas écraser
    // les données avant qu'elles ne soient déplacées
    for (slong i = n - 1; i > 0; i--) {
      error = gr_set(gr_mat_entry_ptr(Res, i, j, ctx), gr_mat_entry_srcptr(M, i - 1, j, ctx), ctx);
      if (error != 0) { return error; }
    }
    // La ligne 0 devient 0
    error = gr_zero(gr_mat_entry_ptr(Res, 0, j, ctx), ctx);
    if (error != 0) { return error; }
  }
  return GR_SUCCESS;
}

int gr_mat_apply_Zt(gr_mat_t Res, gr_mat_t M, gr_ctx_t ctx) {
  int error;
  slong n = gr_mat_nrows(M, ctx);
  slong m = gr_mat_ncols(M, ctx);

  for (slong j = 0; j < m; j++) {
    // Pour Zt (vers le haut), on part du haut vers le bas
    for (slong i = 0; i < n - 1; i++) {
      error = gr_set(gr_mat_entry_ptr(Res, i, j, ctx), gr_mat_entry_srcptr(M, i + 1, j, ctx), ctx);
      if (error != 0) { return error; }
    }
    // La dernière ligne devient 0
    error = gr_zero(gr_mat_entry_ptr(Res, n - 1, j, ctx), ctx);
    if (error != 0) { return error; }
  }
  return GR_SUCCESS;
}
int gr_mat_mul_generator(gr_mat_t G_c, gr_mat_t H_c, gr_mat_t G_a, gr_mat_t H_a, gr_mat_t G_b, gr_mat_t H_b,
                         gr_ctx_t ctx) {
  slong n = gr_mat_nrows(G_a, ctx);
  gr_mat_t W, V, a, b, Tmp, LastCol;
  int error;

  gr_mat_init(W, n, gr_mat_ncols(G_b, ctx), ctx);
  gr_mat_init(V, n, gr_mat_ncols(H_a, ctx), ctx);
  gr_mat_init(a, n, 1, ctx);
  gr_mat_init(b, n, 1, ctx);
  gr_mat_init(Tmp, n, gr_mat_ncols(G_b, ctx), ctx);
  gr_mat_init(LastCol, n, 1, ctx);

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

  // 3. Extraction des colonnes de correction (ZA * e_{n-1} et ZB^T * e_{n-1})
  error = gr_one(gr_mat_entry_ptr(LastCol, n - 1, 0, ctx), ctx);
  if (error != 0) { return error; }

  error = gr_mat_mul_vector(a, G_a, H_a, LastCol, ctx);
  if (error != 0) { return error; }
  error = gr_mat_apply_Z(a, a, ctx);
  if (error != 0) { return error; }

  error = gr_mat_mul_vector(b, H_b, G_b, LastCol, ctx);
  if (error != 0) { return error; }
  error = gr_mat_apply_Z(b, b, ctx);
  if (error != 0) { return error; }

  // 4. Concaténation finale via gr_mat_concat_horizontal
  gr_mat_t G_res, H_res, G_temp, H_temp;
  gr_mat_init(G_res, n, gr_mat_ncols(G_a, ctx) + gr_mat_ncols(W, ctx) + 1, ctx);
  gr_mat_init(G_temp, n, gr_mat_ncols(G_a, ctx) + gr_mat_ncols(W, ctx), ctx);
  gr_mat_init(H_res, n, gr_mat_ncols(V, ctx) + gr_mat_ncols(H_b, ctx) + 1, ctx);
  gr_mat_init(H_temp, n, gr_mat_ncols(V, ctx) + gr_mat_ncols(H_b, ctx), ctx);

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

  error = gr_mat_init_set(G_c, G_res, ctx);
  if (error != 0) { return error; }
  error = gr_mat_init_set(H_c, H_res, ctx);
  if (error != 0) { return error; }

  gr_mat_clear(W, ctx);
  gr_mat_clear(V, ctx);
  gr_mat_clear(a, ctx);
  gr_mat_clear(b, ctx);
  gr_mat_clear(Tmp, ctx);
  gr_mat_clear(LastCol, ctx);
  gr_mat_clear(G_temp, ctx);
  gr_mat_clear(H_temp, ctx);
  gr_mat_clear(G_res, ctx);
  gr_mat_clear(H_res, ctx);

  return GR_SUCCESS;
}
