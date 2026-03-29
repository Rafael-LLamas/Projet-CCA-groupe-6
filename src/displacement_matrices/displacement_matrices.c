#include <stdio.h>
#include <stdlib.h>

#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"
#include "matrix_aux.h"

int gr_mat_displacement_square_safe(gr_mat_t D, gr_mat_t A, gr_ctx_t ctx) {

  /*
   * This is the safest by book definition of the displacement operation
   * seen on the materials.
   *
   * ∇A = A−ZAZ^T
   *
   * The implementation follows the matrix structures of Flint, creating a Z
   * then doing matrix multiplications and substraction to recieve the desired output.
   *
   * We use this generally to base off tests and such.
   */

  int status = GR_SUCCESS;
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
    status |= gr_one(val, ctx);
  }

  status |= gr_mat_transpose(ZT, Z, ctx); // Z^T

  status |= gr_mat_mul(Temp1, A, ZT, ctx); // Temp1 = A * Z^T

  status |= gr_mat_mul(Temp2, Z, Temp1, ctx); // Temp2 = Z * Temp1  (which is Z * A * Z^T)

  status |= gr_mat_sub(D, A, Temp2, ctx); // D = A - Temp2

  gr_mat_clear(Z, ctx);
  gr_mat_clear(ZT, ctx);
  gr_mat_clear(Temp1, ctx);
  gr_mat_clear(Temp2, ctx);
  return status;
}

int gr_mat_displacement(gr_mat_t D, gr_mat_t A, disp_type_t type, gr_ctx_t ctx) {

  /*
   * Computational approach to the previously seen displacement matrix generation.
   *
   * This function is an pointer arithmetic implementation of the previously
   * implemented gr_mat_displacement_square_safe(), we have 2 main options
   * for generation of two different displacement matrices mentioned as Phi-
   * and Phi+ on materials.
   *
   * For displacment operators:
   * Two types can be passed down to the function with the parameter type,
   * the type definitions are located in header file of this file.
   *
   *          DISP_MINUS  -> Phi-
   *          DISP_PLUS   -> Phi+ = ∇A
   *
   * Both will respectfull give the resulting matrix of A−Z^TAZ or A−ZAZ^T.
   */

  int status = GR_SUCCESS;
  slong n = gr_mat_nrows(A, ctx);
  slong m = gr_mat_ncols(A, ctx);

  if (gr_mat_nrows(D, ctx) != n || gr_mat_ncols(D, ctx) != m) return GR_UNABLE;
  gr_ptr ptr_cur, ptr_op, ptr_dest;

  if (type == DISP_MINUS) {

    // inner part of Phi-
    for (slong i = 0; i < n - 1; i++) {
      for (slong j = 0; j < m - 1; j++) {
        ptr_cur = gr_mat_entry_ptr(A, i, j, ctx);
        ptr_op = gr_mat_entry_ptr(A, i + 1, j + 1, ctx);
        ptr_dest = gr_mat_entry_ptr(D, i, j, ctx);
        status |= gr_sub(ptr_dest, ptr_cur, ptr_op, ctx); // D[i,j] = A[i,j] - A[i+1, j+1]
      }
    }

    // copy the rest (last row & column)
    for (slong i = 0; i < n; i++)
      status |= gr_set(gr_mat_entry_ptr(D, i, m - 1, ctx), gr_mat_entry_ptr(A, i, m - 1, ctx), ctx);
    for (slong j = 1; j < m; j++)
      status |= gr_set(gr_mat_entry_ptr(D, n - 1, j, ctx), gr_mat_entry_ptr(A, n - 1, j, ctx), ctx);

  } else { // DISP_PLUS

    // inner part of Phi+
    for (slong i = n - 1; i > 0; i--) {
      for (slong j = 1; j < m; j++) {
        ptr_cur = gr_mat_entry_ptr(A, i, j, ctx);
        ptr_op = gr_mat_entry_ptr(A, i - 1, j - 1, ctx);
        ptr_dest = gr_mat_entry_ptr(D, i, j, ctx);
        status |= gr_sub(ptr_dest, ptr_cur, ptr_op, ctx); // D[i,j] = A[i,j] - A[i-1, j-1]
      }
    }

    // copy the rest (first row & columns)
    for (slong i = 0; i < n; i++) status |= gr_set(gr_mat_entry_ptr(D, i, 0, ctx), gr_mat_entry_ptr(A, i, 0, ctx), ctx);
    for (slong j = 1; j < m; j++) status |= gr_set(gr_mat_entry_ptr(D, 0, j, ctx), gr_mat_entry_ptr(A, 0, j, ctx), ctx);
  }

  return status;
}

int gr_mat_G_H(gr_mat_t G, gr_mat_t H, gr_mat_t A, disp_type_t type, gr_ctx_t ctx) {

  /*
   * Returning the generators of a matrix A,
   * also can be called the rank optimization part of the toeplitz matrices.
   *
   * The displacement matrix of A (<- D) can be factored into G and H^T.
   * The product of G and H^T will give us the displacement matrix back.
   *
   * The mathematics definitions two types of displacement operators:
   *            Phi+ (DISP_PLUS):   D = A - Z * A * Z^T
   *            Phi- (DISP_MINUS):  D = A - Z^T * A * Z
   *
   * Regardless of which operator generated D, the resulting displacement matrix
   * can be expressed as a sum of outer products: D = G * H^T.
   *
   * Because our goal here is strictly to factorize D into a left part (G)
   * and a right part (H^T), we use an LU decomposition
   *
   * (D = P * L * U) for both types without changing the logic. By setting
   * G = P * L and H^T = U, we get G * H^T = D.
   */

  int status = GR_SUCCESS;
  slong m = gr_mat_nrows(A, ctx);
  slong n = gr_mat_ncols(A, ctx);

  slong *P = flint_malloc(m * sizeof(slong)); // Permutation table

  gr_mat_t D;
  gr_mat_init(D, m, n, ctx);

  status |= gr_mat_displacement(D, A, type, ctx); // D <- displacement matrix

  if (status != GR_SUCCESS) {
    gr_mat_clear(D, ctx);
    flint_free(P);
    return status;
  }

  gr_mat_t LU, L, U;
  slong rank;

  gr_mat_init(LU, m, n, ctx);
  gr_mat_init(L, m, m, ctx);
  gr_mat_init(U, m, n, ctx);

  status |= gr_mat_lu(&rank, P, LU, D, 0, ctx); // LU decomposition
  
  if (rank < 1) {
      flint_free(P);
      gr_mat_clear(D, ctx);
      gr_mat_clear(LU, ctx);
      gr_mat_clear(L, ctx);
      gr_mat_clear(U, ctx);
      return GR_UNABLE;
  }

  gr_mat_init(G, m, rank, ctx);
  gr_mat_init(H, n, rank, ctx);

  status |= gr_mat_lu_detach(L, U, LU, ctx); // detach the LU format to L and U

  // copy the values from L to G with correct permutation
  for (slong i = 0; i < m; i++)
    for (slong j = 0; j < rank; j++)
      status |= gr_set(gr_mat_entry_ptr(G, P[i], j, ctx), gr_mat_entry_srcptr(L, i, j, ctx), ctx);

  // copy the values from U to H
  for (slong i = 0; i < rank; i++)
    for (slong j = 0; j < n; j++)
      status |= gr_set(gr_mat_entry_ptr(H, j, i, ctx), gr_mat_entry_srcptr(U, i, j, ctx), ctx);

  flint_free(P);
  gr_mat_clear(D, ctx);
  gr_mat_clear(LU, ctx);
  gr_mat_clear(L, ctx);
  gr_mat_clear(U, ctx);
  return status;
}

int gr_mat_reconstruct_A(gr_mat_t A, gr_mat_t G, gr_mat_t H, disp_type_t type, gr_ctx_t ctx) {

  /* Reconstructs the original A from its generators
   *
   * This function performs the inverse of the displacement operator
   * It is a computational approach to the multiplication of L and U using
   * diagonal pointer arithmetics
   *
   * The reconstruction depends on the displacement type used to generate
   * G and H:
   *
   * DISP_PLUS | Phi+ (Sigma LU)
   * Toeplitz data flows from the top-left corner, to compute cell A[i,j],
   * following the data backwards by walking diagonally UP and LEFT towards
   * the (0,0). We accumulate G[i-x] * H[j-x] and walk until we hit
   * the top or left wall (max_x = min(i, j)).
   *
   * DISP_MINUS | Phi- (Sigma UL)
   * A = Sum( U(rev(x_k)) * L(rev(y_k)) )
   * Toeplitz data flows from the bottom-right corner. We follow the data
   * by walking diagonally DOWN and RIGHT. Accumulate G[i+x] * H[j+x]
   * walking until we hit the bottom or right wall (max_x = min(n-1-i, m-1-j)).
   */

  int status = GR_SUCCESS;
  slong rank = gr_mat_ncols(G, ctx);

  if (gr_mat_nrows(A, ctx) != gr_mat_nrows(G, ctx) || gr_mat_ncols(A, ctx) != gr_mat_nrows(H, ctx)) {
    return GR_UNABLE;
  }

  gr_ptr sum_res = gr_heap_init(ctx); // final value for A[i,j]
  gr_ptr temp = gr_heap_init(ctx);    // g * h

  if (type == DISP_MINUS) {

    for (slong i = 0; i < gr_mat_ncols(A, ctx); i++) {
      for (slong j = 0; j < gr_mat_nrows(A, ctx); j++) {
        status |= gr_zero(sum_res, ctx);
        for (slong k = 0; k < rank; k++) {
          slong max_x = FLINT_MIN(gr_mat_ncols(A, ctx) - 1 - i, gr_mat_nrows(A, ctx) - 1 - j);
          for (slong x = 0; x <= max_x; x++) {
            status |= gr_mul(temp, gr_mat_entry_ptr(G, i + x, k, ctx), gr_mat_entry_ptr(H, j + x, k, ctx), ctx);
            status |= gr_add(sum_res, sum_res, temp, ctx);
          }
        }
        status |= gr_set(gr_mat_entry_ptr(A, i, j, ctx), sum_res, ctx);
      }
    }

  } else { // DISP_PLUS

    for (slong i = 0; i < gr_mat_ncols(A, ctx); i++) {
      for (slong j = 0; j < gr_mat_nrows(A, ctx); j++) {
        status |= gr_zero(sum_res, ctx);
        for (slong k = 0; k < rank; k++) {
          slong max_x = FLINT_MIN(i, j);
          for (slong x = 0; x <= max_x; x++) {
            status |= gr_mul(temp, gr_mat_entry_ptr(G, i - x, k, ctx), gr_mat_entry_ptr(H, j - x, k, ctx), ctx);
            status |= gr_add(sum_res, sum_res, temp, ctx);
          }
        }
        status |= gr_set(gr_mat_entry_ptr(A, i, j, ctx), sum_res, ctx);
      }
    }
  }

  gr_heap_clear(sum_res, ctx);
  gr_heap_clear(temp, ctx);
  return status;
}
