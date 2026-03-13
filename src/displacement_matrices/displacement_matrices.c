#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "matrix_aux.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int gr_mat_displacement_square_safe(gr_mat_t D, gr_mat_t A, gr_ctx_t ctx) {
  slong n = gr_mat_nrows(A, ctx);
  gr_mat_t Z, ZT, Temp1, Temp2;
  gr_mat_init(Z, n, n, ctx);
  gr_mat_init(ZT, n, n, ctx);
  gr_mat_init(Temp1, n, n, ctx);
  gr_mat_init(Temp2, n, n, ctx);
  FLINT_CHECK(gr_mat_zero(Z, ctx));
  gr_ptr val;
  for (slong i = 1; i < n; i++) { // create Z
    val = gr_mat_entry_ptr(Z, i, i - 1, ctx);
    FLINT_CHECK(gr_one(val, ctx));
  }
  FLINT_CHECK(gr_mat_transpose(ZT, Z, ctx));     // Z^T
  FLINT_CHECK(gr_mat_mul(Temp1, A, ZT, ctx));    // Temp1 = A * Z^T
  FLINT_CHECK(gr_mat_mul(Temp2, Z, Temp1, ctx)); // Temp2 = Z * Temp1  (which is Z * A * Z^T)
  FLINT_CHECK(gr_mat_sub(D, A, Temp2, ctx));     // D = A - Temp2
  gr_mat_clear(Z, ctx);
  gr_mat_clear(ZT, ctx);
  gr_mat_clear(Temp1, ctx);
  gr_mat_clear(Temp2, ctx);
  return GR_SUCCESS;
}

int gr_mat_displacement(gr_mat_t D, gr_mat_t A, gr_ctx_t ctx, disp_type_t type) {
  slong n = gr_mat_nrows(A, ctx);
  slong m = gr_mat_ncols(A, ctx);

  if (gr_mat_nrows(D, ctx) != n || gr_mat_ncols(D, ctx) != m) return GR_UNABLE;
  gr_ptr ptr_cur, ptr_op, ptr_dest;

  if (type == DISP_MINUS) {

    for (slong i = 0; i < n - 1; i++) {
      for (slong j = 0; j < m - 1; j++) {
        ptr_cur = gr_mat_entry_ptr(A, i, j, ctx);
        ptr_op = gr_mat_entry_ptr(A, i + 1, j + 1, ctx);
        ptr_dest = gr_mat_entry_ptr(D, i, j, ctx);
        FLINT_CHECK(gr_sub(ptr_dest, ptr_cur, ptr_op, ctx)); // D[i,j] = A[i,j] - A[i+1, j+1]
      }
    }

    // copy the rest (left row & column)
    for (slong i = 0; i < n; i++)
      FLINT_CHECK(gr_set(gr_mat_entry_ptr(D, i, m-1, ctx), gr_mat_entry_ptr(A, i, m-1, ctx), ctx));
    for (slong j = 1; j < m; j++)
      FLINT_CHECK(gr_set(gr_mat_entry_ptr(D, n-1, j, ctx), gr_mat_entry_ptr(A, n-1, j, ctx), ctx));

  } else {

    // inner part
    for (slong i = n - 1; i > 0; i--) {
      for (slong j = 1; j < m; j++) {
        ptr_cur = gr_mat_entry_ptr(A, i, j, ctx);
        ptr_op = gr_mat_entry_ptr(A, i - 1, j - 1, ctx);
        ptr_dest = gr_mat_entry_ptr(D, i, j, ctx);
        FLINT_CHECK(gr_sub(ptr_dest, ptr_cur, ptr_op, ctx)); // D[i,j] = A[i,j] - A[i-1, j-1]
      }
    }

    // copy the rest (first row & columns)
    for (slong i = 0; i < n; i++)
      FLINT_CHECK(gr_set(gr_mat_entry_ptr(D, i, 0, ctx), gr_mat_entry_ptr(A, i, 0, ctx), ctx));
    for (slong j = 1; j < m; j++)
      FLINT_CHECK(gr_set(gr_mat_entry_ptr(D, 0, j, ctx), gr_mat_entry_ptr(A, 0, j, ctx), ctx));
  }

  return GR_SUCCESS;
}

int gr_mat_G_H(gr_mat_t G, gr_mat_t H, gr_mat_t A, gr_ctx_t ctx) {
  slong m = gr_mat_nrows(A, ctx);
  slong n = gr_mat_ncols(A, ctx);
  slong *P = flint_malloc(m * sizeof(slong));
  gr_mat_t D;
  gr_mat_init(D, m, n, ctx);
  int res = gr_mat_displacement(D, A, ctx, DISP_PLUS);
  if (res != GR_SUCCESS) {
    gr_mat_clear(D, ctx);
    flint_free(P);
    return res;
  }
  gr_mat_t LU, L, U;
  slong rank;
  gr_mat_init(LU, m, n, ctx);
  gr_mat_init(L, m, m, ctx);
  gr_mat_init(U, m, n, ctx);
  FLINT_CHECK(gr_mat_lu(&rank, P, LU, D, 0, ctx));
  gr_mat_init(G, m, rank, ctx);
  gr_mat_init(H, n, rank, ctx);
  FLINT_CHECK(gr_mat_lu_detach(L, U, LU, ctx));
  for (slong i = 0; i < m; i++)
    for (slong j = 0; j < rank; j++)
      FLINT_CHECK(gr_set(gr_mat_entry_ptr(G, P[i], j, ctx), gr_mat_entry_srcptr(L, i, j, ctx), ctx));
  for (slong i = 0; i < rank; i++)
    for (slong j = 0; j < n; j++)
      FLINT_CHECK(gr_set(gr_mat_entry_ptr(H, j, i, ctx), gr_mat_entry_srcptr(U, i, j, ctx), ctx));
  flint_free(P);
  gr_mat_clear(D, ctx);
  gr_mat_clear(LU, ctx);
  gr_mat_clear(L, ctx);
  gr_mat_clear(U, ctx);
  return GR_SUCCESS;
}

int gr_mat_reconstruct_A_safe(gr_mat_t A, gr_mat_t G, gr_mat_t H, gr_ctx_t ctx) {
  slong rank = gr_mat_ncols(G, ctx);
  if (gr_mat_nrows(A, ctx) != gr_mat_nrows(G, ctx) || gr_mat_ncols(A, ctx) != gr_mat_nrows(H, ctx)) {
    return GR_UNABLE;
  }
  gr_ptr sum_res = gr_heap_init(ctx);                // final value for A[i,j]
  gr_ptr temp = gr_heap_init(ctx);                   // g * h
  for (slong i = 0; i < gr_mat_nrows(G, ctx); i++) { // for every element of A
    for (slong j = 0; j < gr_mat_nrows(H, ctx); j++) {
      FLINT_CHECK(gr_zero(sum_res, ctx));
      for (slong k = 0; k < rank; k++) { // Sigma - calculate this directly (L_k * U_k)[i,j]
        slong minij = FLINT_MIN(i, j);
        for (slong x = 0; x <= minij; x++) {
          // L[i, x] = G[i-x, k]
          // U[x, j] = H[j-x, k]
          FLINT_CHECK(gr_mul(temp, gr_mat_entry_ptr(G, i - x, k, ctx), gr_mat_entry_ptr(H, j - x, k, ctx), ctx));
          FLINT_CHECK(gr_add(sum_res, sum_res, temp, ctx));
        }
      }
      FLINT_CHECK(gr_set(gr_mat_entry_ptr(A, i, j, ctx), sum_res, ctx));
    }
  }
  gr_heap_clear(sum_res, ctx);
  gr_heap_clear(temp, ctx);
  return GR_SUCCESS;
}
