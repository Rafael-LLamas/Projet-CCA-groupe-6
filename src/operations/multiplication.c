#include "matrix_aux.h"

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include <stdlib.h>
#include <time.h>

int gr_mat_mul_vector(gr_mat_t Res, gr_mat_t G, gr_mat_t H, gr_mat_t X, gr_ctx_t ctx) {
  int error = GR_SUCCESS;
  slong n = gr_mat_nrows(G, ctx), m = gr_mat_nrows(H, ctx);
  slong alpha = gr_mat_ncols(G, ctx), m_cols = gr_mat_ncols(X, ctx);

  gr_poly_t pg, ph, px, p_tmp;
  gr_poly_init(pg, ctx);
  gr_poly_init(ph, ctx);
  gr_poly_init(px, ctx);
  gr_poly_init(p_tmp, ctx);

  error = gr_mat_zero(Res, ctx);
  if (error) {
    gr_poly_clear(pg, ctx);
    gr_poly_clear(ph, ctx);
    gr_poly_clear(px, ctx);
    gr_poly_clear(p_tmp, ctx);
  }

  for (slong j = 0; j < m_cols; j++) {
    gr_poly_fit_length(px, m, ctx);
    for (slong r = 0; r < m; r++) {
      error = gr_set(gr_poly_coeff_ptr(px, r, ctx), gr_mat_entry_srcptr(X, r, j, ctx), ctx);
      if (error) {
        gr_poly_clear(pg, ctx);
        gr_poly_clear(ph, ctx);
        gr_poly_clear(px, ctx);
        gr_poly_clear(p_tmp, ctx);
      }
    }
    _gr_poly_set_length(px, m, ctx);

    for (slong k = 0; k < alpha; k++) {
      gr_poly_fit_length(pg, n, ctx);
      for (slong r = 0; r < n; r++) {
        error = gr_set(gr_poly_coeff_ptr(pg, r, ctx), gr_mat_entry_srcptr(G, r, k, ctx), ctx);
        if (error) {
          gr_poly_clear(pg, ctx);
          gr_poly_clear(ph, ctx);
          gr_poly_clear(px, ctx);
          gr_poly_clear(p_tmp, ctx);
        }
      }
      _gr_poly_set_length(pg, n, ctx);

      gr_poly_fit_length(ph, m, ctx);
      for (slong r = 0; r < m; r++) {
        error = gr_set(gr_poly_coeff_ptr(ph, r, ctx), gr_mat_entry_srcptr(H, r, k, ctx), ctx);
        if (error) {
          gr_poly_clear(pg, ctx);
          gr_poly_clear(ph, ctx);
          gr_poly_clear(px, ctx);
          gr_poly_clear(p_tmp, ctx);
        }
      }
      _gr_poly_set_length(ph, m, ctx);

      error = gr_poly_reverse(ph, ph, m, ctx);
      if (error) {
        gr_poly_clear(pg, ctx);
        gr_poly_clear(ph, ctx);
        gr_poly_clear(px, ctx);
        gr_poly_clear(p_tmp, ctx);
      }
      error = gr_poly_mul(p_tmp, ph, px, ctx);
      if (error) {
        gr_poly_clear(pg, ctx);
        gr_poly_clear(ph, ctx);
        gr_poly_clear(px, ctx);
        gr_poly_clear(p_tmp, ctx);
      }

      error = gr_poly_shift_right(p_tmp, p_tmp, m - 1, ctx);
      if (error) {
        gr_poly_clear(pg, ctx);
        gr_poly_clear(ph, ctx);
        gr_poly_clear(px, ctx);
        gr_poly_clear(p_tmp, ctx);
      }
      _gr_poly_set_length(p_tmp, FLINT_MIN(gr_poly_length(p_tmp, ctx), m), ctx);

      error = gr_poly_mul(p_tmp, pg, p_tmp, ctx);
      if (error) {
        gr_poly_clear(pg, ctx);
        gr_poly_clear(ph, ctx);
        gr_poly_clear(px, ctx);
        gr_poly_clear(p_tmp, ctx);
      }

      if (gr_poly_length(p_tmp, ctx) > n) _gr_poly_set_length(p_tmp, n, ctx);

      for (slong i = 0; i < gr_poly_length(p_tmp, ctx); i++) {
        error = gr_add(gr_mat_entry_ptr(Res, i, j, ctx), gr_mat_entry_ptr(Res, i, j, ctx),
                       gr_poly_coeff_srcptr(p_tmp, i, ctx), ctx);
        if (error) {
          gr_poly_clear(pg, ctx);
          gr_poly_clear(ph, ctx);
          gr_poly_clear(px, ctx);
          gr_poly_clear(p_tmp, ctx);
        }
      }
    }
  }

  gr_poly_clear(pg, ctx);
  gr_poly_clear(ph, ctx);
  gr_poly_clear(px, ctx);
  gr_poly_clear(p_tmp, ctx);
  return error;
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

  // 1. W = Z * A * Z^T * G_b
  error = gr_mat_apply_Zt(Tmp, G_b, ctx);

  if (error != 0) {
    gr_mat_clear(W, ctx);
    gr_mat_clear(V, ctx);
    gr_mat_clear(a, ctx);
    gr_mat_clear(b, ctx);
    gr_mat_clear(Tmp, ctx);

    return error;
  }
  error = gr_mat_mul_vector(W, G_a, H_a, Tmp, ctx);

  if (error != 0) {
    gr_mat_clear(W, ctx);
    gr_mat_clear(V, ctx);
    gr_mat_clear(a, ctx);
    gr_mat_clear(b, ctx);
    gr_mat_clear(Tmp, ctx);

    return error;
  }
  error = gr_mat_apply_Z(W, W, ctx);

  if (error != 0) {
    gr_mat_clear(W, ctx);
    gr_mat_clear(V, ctx);
    gr_mat_clear(a, ctx);
    gr_mat_clear(b, ctx);
    gr_mat_clear(Tmp, ctx);

    return error;
  }

  // 2. V = B^T * H_a (Le générateur de B^T est {H_b, G_b})
  error = gr_mat_mul_vector(V, H_b, G_b, H_a, ctx);

  if (error != 0) {
    gr_mat_clear(W, ctx);
    gr_mat_clear(V, ctx);
    gr_mat_clear(a, ctx);
    gr_mat_clear(b, ctx);
    gr_mat_clear(Tmp, ctx);

    return error;
  }

  // 3. Correction colonnes
  // Vecteur e_{m-1} pour la dimension de A
  // LastCol_m doit être e_{nrows(H_a) - 1}
  slong m_a = gr_mat_nrows(H_a, ctx);
  gr_mat_init(LastCol_m, m_a, 1, ctx);
  error = gr_mat_zero(LastCol_m, ctx);
  if (error != 0) {
    gr_mat_clear(W, ctx);
    gr_mat_clear(V, ctx);
    gr_mat_clear(a, ctx);
    gr_mat_clear(b, ctx);
    gr_mat_clear(Tmp, ctx);
    gr_mat_clear(LastCol_m, ctx);

    return error;
  }

  error = gr_one(gr_mat_entry_ptr(LastCol_m, m_a - 1, 0, ctx), ctx);

  if (error != 0) {
    gr_mat_clear(W, ctx);
    gr_mat_clear(V, ctx);
    gr_mat_clear(a, ctx);
    gr_mat_clear(b, ctx);
    gr_mat_clear(Tmp, ctx);
    gr_mat_clear(LastCol_m, ctx);

    return error;
  }

  error = gr_mat_mul_vector(a, G_a, H_a, LastCol_m, ctx);
  if (error != 0) {
    gr_mat_clear(W, ctx);
    gr_mat_clear(V, ctx);
    gr_mat_clear(a, ctx);
    gr_mat_clear(b, ctx);
    gr_mat_clear(Tmp, ctx);
    gr_mat_clear(LastCol_m, ctx);

    return error;
  }
  error = gr_mat_apply_Z(a, a, ctx);

  if (error != 0) {
    gr_mat_clear(W, ctx);
    gr_mat_clear(V, ctx);
    gr_mat_clear(a, ctx);
    gr_mat_clear(b, ctx);
    gr_mat_clear(Tmp, ctx);
    gr_mat_clear(LastCol_m, ctx);
    return error;
  }

  // Vecteur e_{m-1} pour la dimension de B^T
  // B^T a pour générateurs {H_b, G_b}, donc m = nrows(G_b)
  slong m_bt = gr_mat_nrows(G_b, ctx);
  gr_mat_init(LastCol_k, m_bt, 1, ctx);

  error = gr_mat_zero(LastCol_k, ctx);
  if (error != 0) {
    gr_mat_clear(W, ctx);
    gr_mat_clear(V, ctx);
    gr_mat_clear(a, ctx);
    gr_mat_clear(b, ctx);
    gr_mat_clear(Tmp, ctx);
    gr_mat_clear(LastCol_m, ctx);
    gr_mat_clear(LastCol_k, ctx);
    return error;
  }

  error = gr_one(gr_mat_entry_ptr(LastCol_k, m_bt - 1, 0, ctx), ctx);
  if (error != 0) {
    gr_mat_clear(W, ctx);
    gr_mat_clear(V, ctx);
    gr_mat_clear(a, ctx);
    gr_mat_clear(b, ctx);
    gr_mat_clear(Tmp, ctx);
    gr_mat_clear(LastCol_m, ctx);
    gr_mat_clear(LastCol_k, ctx);
    return error;
  }

  error = gr_mat_mul_vector(b, H_b, G_b, LastCol_k, ctx);

  if (error != 0) {
    gr_mat_clear(W, ctx);
    gr_mat_clear(V, ctx);
    gr_mat_clear(a, ctx);
    gr_mat_clear(b, ctx);
    gr_mat_clear(Tmp, ctx);
    gr_mat_clear(LastCol_m, ctx);
    gr_mat_clear(LastCol_k, ctx);
    return error;
  }
  error = gr_mat_apply_Z(b, b, ctx);

  if (error != 0) {
    gr_mat_clear(W, ctx);
    gr_mat_clear(V, ctx);
    gr_mat_clear(a, ctx);
    gr_mat_clear(b, ctx);
    gr_mat_clear(Tmp, ctx);
    gr_mat_clear(LastCol_m, ctx);
    gr_mat_clear(LastCol_k, ctx);
    return error;
  }

  gr_mat_clear(Tmp, ctx);
  gr_mat_clear(LastCol_m, ctx);
  gr_mat_clear(LastCol_k, ctx);

  // 4. Concaténation (G_c: n x (ra+rb+1), H_c: k x (ra+rb+1))
  gr_mat_t G_temp, H_temp;
  gr_mat_init(G_c, n, gr_mat_ncols(G_a, ctx) + gr_mat_ncols(W, ctx) + 1, ctx);
  gr_mat_init(G_temp, n, gr_mat_ncols(G_a, ctx) + gr_mat_ncols(W, ctx), ctx);
  gr_mat_init(H_c, k, gr_mat_ncols(V, ctx) + gr_mat_ncols(H_b, ctx) + 1, ctx);
  gr_mat_init(H_temp, k, gr_mat_ncols(V, ctx) + gr_mat_ncols(H_b, ctx), ctx);

  error = gr_mat_concat_horizontal(G_temp, G_a, W, ctx);

  if (error != 0) {
    gr_mat_clear(W, ctx);
    gr_mat_clear(V, ctx);
    gr_mat_clear(a, ctx);
    gr_mat_clear(b, ctx);
    gr_mat_clear(G_temp, ctx);
    gr_mat_clear(H_temp, ctx);
    return error;
  }
  error = gr_mat_concat_horizontal(G_c, G_temp, a, ctx);

  if (error != 0) {
    gr_mat_clear(W, ctx);
    gr_mat_clear(V, ctx);
    gr_mat_clear(a, ctx);
    gr_mat_clear(b, ctx);
    gr_mat_clear(G_temp, ctx);
    gr_mat_clear(H_temp, ctx);
    return error;
  }

  error = gr_mat_neg(b, b, ctx);

  if (error != 0) {
    gr_mat_clear(W, ctx);
    gr_mat_clear(V, ctx);
    gr_mat_clear(a, ctx);
    gr_mat_clear(b, ctx);
    gr_mat_clear(G_temp, ctx);
    gr_mat_clear(H_temp, ctx);
    return error;
  }
  error = gr_mat_concat_horizontal(H_temp, V, H_b, ctx);

  if (error != 0) {
    gr_mat_clear(W, ctx);
    gr_mat_clear(V, ctx);
    gr_mat_clear(a, ctx);
    gr_mat_clear(b, ctx);
    gr_mat_clear(G_temp, ctx);
    gr_mat_clear(H_temp, ctx);
    return error;
  }
  error = gr_mat_concat_horizontal(H_c, H_temp, b, ctx);

  if (error != 0) {
    gr_mat_clear(W, ctx);
    gr_mat_clear(V, ctx);
    gr_mat_clear(a, ctx);
    gr_mat_clear(b, ctx);
    gr_mat_clear(G_temp, ctx);
    gr_mat_clear(H_temp, ctx);
    return error;
  }

  gr_mat_clear(W, ctx);
  gr_mat_clear(V, ctx);
  gr_mat_clear(a, ctx);
  gr_mat_clear(b, ctx);
  gr_mat_clear(G_temp, ctx);
  gr_mat_clear(H_temp, ctx);

  return GR_SUCCESS;
}