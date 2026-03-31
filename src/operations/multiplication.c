#include "matrix_aux.h"

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include <stdlib.h>
#include <time.h>

int gr_mat_mul_vector(gr_mat_t Res, gr_mat_t G, gr_mat_t H, gr_mat_t X, gr_ctx_t ctx) {
  int error = GR_SUCCESS;
  slong n = gr_mat_nrows(G, ctx);
  slong m = gr_mat_nrows(H, ctx);
  slong alpha = gr_mat_ncols(G, ctx);
  slong m_cols = gr_mat_ncols(X, ctx);

  error = gr_mat_zero(Res, ctx);
  if (error != GR_SUCCESS) return error;

  gr_ptr tmp = gr_heap_init(ctx);
  gr_ptr a_ij = gr_heap_init(ctx); // valeur de A[i,jp]
  gr_ptr acc = gr_heap_init(ctx);  // accumulateur pour Res[i,j]

  for (slong j = 0; j < m_cols; j++) {
    for (slong i = 0; i < n; i++) {

      error = gr_zero(acc, ctx);
      if (error != GR_SUCCESS) goto cleanup;

      for (slong jp = 0; jp < m; jp++) {

        // Calcule A[i,jp] = sum_k sum_{x=0}^{min(i,jp)} G[i-x,k]*H[jp-x,k]
        error = gr_zero(a_ij, ctx);
        if (error != GR_SUCCESS) goto cleanup;

        slong minijp = FLINT_MIN(i, jp);
        for (slong x = 0; x <= minijp; x++) {
          for (slong k = 0; k < alpha; k++) {
            error = gr_mul(tmp, gr_mat_entry_srcptr(G, i - x, k, ctx), gr_mat_entry_srcptr(H, jp - x, k, ctx), ctx);
            if (error != GR_SUCCESS) goto cleanup;
            error = gr_add(a_ij, a_ij, tmp, ctx);
            if (error != GR_SUCCESS) goto cleanup;
          }
        }

        // acc += A[i,jp] * X[jp,j]
        error = gr_mul(tmp, a_ij, gr_mat_entry_srcptr(X, jp, j, ctx), ctx);
        if (error != GR_SUCCESS) goto cleanup;
        error = gr_add(acc, acc, tmp, ctx);
        if (error != GR_SUCCESS) goto cleanup;
      }

      error = gr_set(gr_mat_entry_ptr(Res, i, j, ctx), acc, ctx);
      if (error != GR_SUCCESS) goto cleanup;
    }
  }

cleanup:
  gr_heap_clear(tmp, ctx);
  gr_heap_clear(a_ij, ctx);
  gr_heap_clear(acc, ctx);
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
  flint_printf("Matrice Tmp = \n");
  gr_mat_print(Tmp, ctx);
  flint_printf("\n");

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
  error = gr_mat_mul_vector(W, G_a, H_a, Tmp, ctx);

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
  error = gr_mat_apply_Z(W, W, ctx);
  flint_printf("Matrice W = \n");
  gr_mat_print(W, ctx);
  flint_printf("\n");

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

  // 2. V = B^T * H_a (Le générateur de B^T est {H_b, G_b})
  error = gr_mat_mul_vector(V, H_b, G_b, H_a, ctx);
  flint_printf("Matrice v = \n");
  gr_mat_print(V, ctx);
  flint_printf("\n");

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

  // 3. Correction colonnes
  // Vecteur e_{m-1} pour la dimension de A
  // LastCol_m doit être e_{nrows(H_a) - 1}
  slong m_a = gr_mat_nrows(H_a, ctx);
  gr_mat_init(LastCol_m, m_a, 1, ctx);
  gr_mat_zero(LastCol_m, ctx);
  gr_one(gr_mat_entry_ptr(LastCol_m, m_a - 1, 0, ctx), ctx);
  flint_printf("Matrice LastCol_m = \n");
  gr_mat_print(LastCol_m, ctx);
  flint_printf("\n");
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

  error = gr_mat_mul_vector(a, G_a, H_a, LastCol_m, ctx);
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
  error = gr_mat_apply_Z(a, a, ctx);
  flint_printf("Matrice a = \n");
  gr_mat_print(a, ctx);
  flint_printf("\n");

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

  // Vecteur e_{m-1} pour la dimension de B^T
  // B^T a pour générateurs {H_b, G_b}, donc m = nrows(G_b)
  slong m_bt = gr_mat_nrows(G_b, ctx);
  gr_mat_init(LastCol_k, m_bt, 1, ctx); // ← taille m_bt, pas k

  error = gr_mat_zero(LastCol_k, ctx);
  if (error != GR_SUCCESS) { /* cleanup */
    return error;
  }

  error = gr_one(gr_mat_entry_ptr(LastCol_k, m_bt - 1, 0, ctx), ctx);
  if (error != GR_SUCCESS) { /* cleanup */
    return error;
  }
  flint_printf("Matrice LastCol_k = \n");
  gr_mat_print(LastCol_k, ctx);
  flint_printf("\n");
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
  flint_printf("Matrice b = \n");
  gr_mat_print(b, ctx);
  flint_printf("\n");

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
  flint_printf("Matrice G_temp = \n");
  gr_mat_print(G_temp, ctx);
  flint_printf("\n");
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
  flint_printf("Matrice G_c = \n");
  gr_mat_print(G_c, ctx);
  flint_printf("\n");

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
  flint_printf("Matrice neg b = \n");
  gr_mat_print(b, ctx);
  flint_printf("\n");
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
  flint_printf("Matrice H_temp = \n");
  gr_mat_print(H_temp, ctx);
  flint_printf("\n");
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
  flint_printf("Matrice H_c = \n");
  gr_mat_print(H_c, ctx);
  flint_printf("\n");

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