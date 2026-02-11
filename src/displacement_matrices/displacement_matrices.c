#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "matrix_aux.h"
#include "random_toeplitz.h"

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

int gr_mat_displacement(gr_mat_t D, gr_mat_t A, gr_ctx_t ctx) {
  slong n = gr_mat_nrows(A, ctx);
  slong m = gr_mat_ncols(A, ctx);
  if (gr_mat_nrows(D, ctx) != n || gr_mat_ncols(D, ctx) != m) return GR_UNABLE;
  gr_ptr ptr_cur, ptr_prev, ptr_dest;

  // inner part
  for (slong i = n - 1; i > 0; i--) {
    for (slong j = 1; j < m; j++) {
      ptr_cur = gr_mat_entry_ptr(A, i, j, ctx);
      ptr_prev = gr_mat_entry_ptr(A, i - 1, j - 1, ctx);
      ptr_dest = gr_mat_entry_ptr(D, i, j, ctx);
      FLINT_CHECK(gr_sub(ptr_dest, ptr_cur, ptr_prev, ctx)); // D[i,j] = A[i,j] - A[i-1, j-1]
    }
  }
  // copy the rest (first row & columns)
  for (slong i = 0; i < n; i++)
    FLINT_CHECK(gr_set(gr_mat_entry_ptr(D, i, 0, ctx), gr_mat_entry_ptr(A, i, 0, ctx), ctx));
  for (slong j = 1; j < m; j++)
    FLINT_CHECK(gr_set(gr_mat_entry_ptr(D, 0, j, ctx), gr_mat_entry_ptr(A, 0, j, ctx), ctx));
  return GR_SUCCESS;
}

int gr_mat_G_H(gr_mat_t G, gr_mat_t H, gr_mat_t A, slong *rank, gr_ctx_t ctx) {
  slong m = gr_mat_nrows(A, ctx);
  slong n = gr_mat_ncols(A, ctx);
  slong rank_displacement;
  slong *P = flint_malloc(m * sizeof(slong));
  gr_mat_displacement(A, A, ctx);
  FLINT_CHECK(gr_mat_rank(&rank_displacement, A, ctx));
  if (rank_displacement != *rank) {
    if (*rank == -1) {
      gr_mat_clear(H, ctx);
      gr_mat_clear(G, ctx);
    }
    gr_mat_init(G, m, rank_displacement, ctx);
    gr_mat_init(H, n, rank_displacement, ctx);
    *rank = rank_displacement;
  }
  gr_mat_t LU, L, U;
  gr_mat_init(LU, m, n, ctx);
  gr_mat_init(L, m, m, ctx);
  gr_mat_init(U, m, n, ctx);
  FLINT_CHECK(gr_mat_lu(rank, P, LU, A, 0, ctx));
  FLINT_CHECK(gr_mat_lu_detach(L, U, LU, ctx));
  for (slong i = 0; i < m; i++) // extract G
    for (slong j = 0; j < *rank; j++)
      FLINT_CHECK(gr_set(gr_mat_entry_ptr(G, i, j, ctx), gr_mat_entry_srcptr(L, i, j, ctx), ctx));
  for (slong i = 0; i < *rank; i++) // extract H
    for (slong j = 0; j < n; j++)
      FLINT_CHECK(gr_set(gr_mat_entry_ptr(H, j, i, ctx), gr_mat_entry_srcptr(U, i, j, ctx), ctx));
  // for (slong i = 0; i < gr_mat_nrows(G, ctx); i++) FLINT_CHECK(gr_mat_move_row(G, i, P[i], ctx));
  flint_free(P);
  gr_mat_clear(LU, ctx);
  gr_mat_clear(L, ctx);
  gr_mat_clear(U, ctx);
  return GR_SUCCESS;
}

int gr_mat_reconstruct_A_safe(gr_mat_t *A, gr_mat_t G, gr_mat_t H, gr_ctx_t ctx) {
  slong m = gr_mat_nrows(G, ctx);
  slong n = gr_mat_nrows(H, ctx);
  slong rank = gr_mat_ncols(G, ctx);
  gr_mat_init(*A, m, n, ctx); // need to discuss this later
  gr_ptr sum_res = gr_heap_init(ctx); // final value for A[i,j]
  gr_ptr temp = gr_heap_init(ctx);    // g * h
  for (slong i = 0; i < m; i++) {     // for every element of A
    for (slong j = 0; j < n; j++) {
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
      FLINT_CHECK(gr_set(gr_mat_entry_ptr(*A, i, j, ctx), sum_res, ctx));
    }
  }
  gr_heap_clear(sum_res, ctx);
  gr_heap_clear(temp, ctx);
  return GR_SUCCESS;
}

int test_displacement_matrices() {
  flint_printf("*----------* Displacement Matrix Test *----------*\n");
  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, GNMOD);
  int res = GR_SUCCESS;

  flint_printf(":-------: Manual nxn Toeplitz Matrix Deplacement Test :-------:\n");
  {
    gr_mat_t A, D;
    gr_mat_init(A, 3, 3, ctx);
    gr_mat_init(D, 3, 3, ctx);
    // feast your eyes
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 0, ctx), 1, ctx)); // Row 0: 1, 3, 4
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 1, ctx), 3, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 2, ctx), 4, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 0, ctx), 8, ctx)); // Row 1: 8, 1, 3
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 1, ctx), 1, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 2, ctx), 3, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 0, ctx), 9, ctx)); // Row 2: 9, 8, 1
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 1, ctx), 8, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 2, ctx), 1, ctx));
    flint_printf("Input Manual Matrix A:\n");
    gr_mat_print(A, ctx);
    flint_printf("\n");
    flint_printf("Input LU Matrix A:\n");
    gr_mat_displacement(D, A, ctx);
    flint_printf("Displacement Matrix (∇A):\n");
    gr_mat_print(D, ctx);
    flint_printf("\n");
    gr_mat_clear(A, ctx);
    gr_mat_clear(D, ctx);
  }

  flint_printf("\n:-------: Manual nxm Matrix Test :-------:\n");
  {
    gr_mat_t A, D;
    gr_mat_init(A, 4, 3, ctx);
    gr_mat_init(D, 4, 3, ctx);
    // (i am too tired to do a double triangular loop future me if ur reading this... dew it -palpatine)
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 0, ctx), 1, ctx)); // Row 0: 1, 2, 3
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 1, ctx), 2, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 2, ctx), 3, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 0, ctx), 4, ctx)); // Row 1: 4, 1, 2
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 1, ctx), 1, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 2, ctx), 2, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 0, ctx), 5, ctx)); // Row 2: 5, 4, 9
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 1, ctx), 4, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 2, ctx), 9, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 3, 0, ctx), 6, ctx)); // Row 3: 6, 5, 4
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 3, 1, ctx), 5, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 3, 2, ctx), 4, ctx));
    flint_printf("Input Manual Matrix A:\n");
    gr_mat_print(A, ctx);
    flint_printf("\n");
    gr_mat_displacement(D, A, ctx);
    flint_printf("Displacement Matrix (∇A):\n");
    gr_mat_print(D, ctx);
    flint_printf("\n");
    slong rank;
    slong *P = flint_malloc(4 * sizeof(slong));
    gr_mat_t LU;
    gr_mat_init(LU, 4, 3, ctx);
    if (gr_mat_lu(&rank, P, LU, D, 0, ctx) == GR_SUCCESS) {
      flint_printf("LU Decomposition successful.\n");
      gr_mat_print(LU, ctx);
    } else
      flint_printf("LU Decomposition failed.\n");
    flint_free(P);
    gr_mat_clear(LU, ctx);
    gr_mat_clear(A, ctx);
    gr_mat_clear(D, ctx);
  }

  flint_printf("\n\n:-------: Random Large nxn Matrix Execution Time Test :-------:\n");
  {
    slong N = 500;
    flint_printf("Generating %d x %d matrix...\n", N, N);
    gr_mat_t A, D1, D2;
    gr_mat_init(A, N, N, ctx);
    gr_mat_init(D1, N, N, ctx);
    gr_mat_init(D2, N, N, ctx);
    flint_rand_t rand_state;
    flint_rand_init(rand_state);
    FLINT_CHECK(gr_mat_randtest(A, rand_state, ctx));
    // flint_printf("Random Matrix A:\n");
    // gr_mat_print(A, ctx);
    clock_t start, end;
    double time_slow, time_fast;
    flint_printf("\n\nRunning Safe Method (Matrix Mul)...\n"); // matrix multiplication
    start = clock();
    gr_mat_displacement_square_safe(D1, A, ctx);
    // flint_printf("Safe Book result D1:\n");
    // gr_mat_print(D1, ctx);
    end = clock();
    time_slow = ((double)(end - start)) / CLOCKS_PER_SEC;
    flint_printf("\nUnoptimized Time: %f seconds\n", time_slow);
    flint_printf("Running Optimized Method (Element-wise)...\n"); // faster one
    start = clock();
    gr_mat_displacement(D2, A, ctx);
    // flint_printf("Safe Book result D2:\n");
    // gr_mat_print(D2, ctx);
    end = clock();
    time_fast = ((double)(end - start)) / CLOCKS_PER_SEC;
    flint_printf("\nOptimized Time: %f seconds\n", time_fast);
    if (gr_mat_equal(D1, D2, ctx) == T_TRUE)
      flint_printf("[SUCCESS] Both results are the same\n");
    else {
      flint_printf("[FAILED] Results are different!\n");
      res = GR_TEST_FAIL;
    }
    gr_mat_clear(A, ctx);
    gr_mat_clear(D1, ctx);
    gr_mat_clear(D2, ctx);
  }

  flint_printf("\n\n:-------: A Toeplitz matrix -> G & H Test :-------:\n");
  {
    gr_mat_t A, B, T, G, H, HT;
    slong m = 7, n = 7, rank = 0;
    slong *P = flint_malloc(sizeof(slong) * m);
    flint_rand_t state;
    gr_ctx_init_nmod(ctx, GNMOD);
    flint_rand_init(state);
    gr_mat_init(A, m, n, ctx);
    gr_mat_init(B, m, n, ctx);
    gr_mat_init(T, m, n, ctx);
    random_toeplitz(A, n, m, state, ctx);
    FLINT_CHECK(gr_mat_set(T, A, ctx));
    flint_printf("Original matrix A:\n");
    gr_mat_print(A, ctx);
    gr_mat_G_H(G, H, A, &rank, ctx);
    flint_printf("\nGenerator G (%ldx%ld):\n", m, rank);
    gr_mat_print(G, ctx);
    flint_printf("\nGenerator H (%ldx%ld):\n", n, rank);
    gr_mat_print(H, ctx);
    gr_mat_init(HT, rank, n, ctx); // reconstruct B = G * H^T
    FLINT_CHECK(gr_mat_transpose(HT, H, ctx));
    FLINT_CHECK(gr_mat_mul(B, G, HT, ctx));
    flint_printf("\nReconstructed Matrix B:\n");
    gr_mat_print(B, ctx);
    flint_printf("\nShould Be equal to displacement matrix of A:\n");
    FLINT_CHECK(gr_mat_displacement(T, T, ctx));
    gr_mat_print(T, ctx);
    flint_free(P);
    gr_mat_clear(A, ctx);
    gr_mat_clear(T, ctx);
    gr_mat_clear(G, ctx);
    gr_mat_clear(H, ctx);
    gr_mat_clear(HT, ctx);
    gr_mat_clear(B, ctx);
  }

  flint_printf("\n\n:-------: A Quasi Toeplitz matrix -> G & H Test :-------:\n");
  {
    gr_mat_t A, B, T, G, H, HT;
    slong m = 7, n = 7, rank = 0;
    slong *P = flint_malloc(sizeof(slong) * m);
    flint_rand_t state;
    gr_ctx_init_nmod(ctx, GNMOD);
    flint_rand_init(state);
    gr_mat_init(A, m, n, ctx);
    gr_mat_init(B, m, n, ctx);
    gr_mat_init(T, m, n, ctx);
    random_quasi_toeplitz(A, n, m, state, ctx);
    FLINT_CHECK(gr_mat_set(T, A, ctx));
    flint_printf("Original matrix A:\n");
    gr_mat_print(A, ctx);
    gr_mat_G_H(G, H, A, &rank, ctx);
    flint_printf("\nGenerator G (%ldx%ld):\n", m, rank);
    gr_mat_print(G, ctx);
    flint_printf("\nGenerator H (%ldx%ld):\n", n, rank);
    gr_mat_print(H, ctx);
    gr_mat_init(HT, rank, n, ctx); // reconstruct B = G * H^T
    FLINT_CHECK(gr_mat_transpose(HT, H, ctx));
    FLINT_CHECK(gr_mat_mul(B, G, HT, ctx));
    flint_printf("\nReconstructed Matrix B:\n");
    gr_mat_print(B, ctx);
    flint_printf("\nShould Be equal to displacement matrix of A:\n");
    FLINT_CHECK(gr_mat_displacement(T, T, ctx));
    gr_mat_print(T, ctx);
    flint_free(P);
    gr_mat_clear(A, ctx);
    gr_mat_clear(T, ctx);
    gr_mat_clear(G, ctx);
    gr_mat_clear(H, ctx);
    gr_mat_clear(HT, ctx);
    gr_mat_clear(B, ctx);
  }

  flint_printf("\n\n:-------: Reconstruction Test :-------:\n");
  {
    gr_mat_t A, A_rec, G, H;
    slong m = 7, n = 7, rank = 0;
    flint_rand_t state;
    gr_ctx_init_nmod(ctx, GNMOD);
    flint_rand_init(state);
    gr_mat_init(A, m, n, ctx);
    random_quasi_toeplitz(A, n, m, state, ctx);
    flint_printf("Original Matrix A:\n");
    gr_mat_print(A, ctx);
    gr_mat_G_H(G, H, A, &rank, ctx);
    FLINT_CHECK(gr_mat_reconstruct_A_safe(&A_rec, G, H, ctx));
    flint_printf("\nReconstructed Matrix:\n");
    gr_mat_print(A_rec, ctx);
    gr_mat_clear(A_rec, ctx);
    gr_mat_clear(A, ctx);
    gr_mat_clear(G, ctx);
    gr_mat_clear(H, ctx);
    flint_rand_clear(state);
  }
  gr_ctx_clear(ctx);
  return res;
}