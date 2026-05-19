#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"
#include "flint/ulong_extras.h"
#include "matrix_aux.h"
#include "random_toeplitz.h"

int test_toeplitz_deplacement() {
  int i = 0;
  flint_rand_t state;
  int error = GR_SUCCESS;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  while (i < 50) {
    slong n = n_randint(state, 100);
    slong m = n_randint(state, 100);
    gr_ctx_t ctx;
    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
    gr_mat_t A, D;
    gr_mat_init(A, n, m, ctx);
    gr_mat_init(D, n, m, ctx);
    error = gr_mat_random_toeplitz(A, state, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to create a Toeplitz A, FIX IT FAST\n");
      gr_mat_clear(A, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_displacement(D, A, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to create the Displacement\n");
      gr_mat_clear(A, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    for (int i = 1; i < gr_mat_nrows(A, ctx); i++) {
      for (int j = 1; j < gr_mat_ncols(A, ctx); j++) {
        error = gr_set_ui(gr_mat_entry_ptr(A, i, j, ctx), 0, ctx);
        if (error != GR_SUCCESS) {
          flint_printf("Failed to put 0\n");
          gr_mat_clear(A, ctx);
          gr_mat_clear(D, ctx);
          gr_ctx_clear(ctx);
          flint_rand_clear(state);
          return error;
        }
      }
    }

    if (gr_mat_equal(A, D, ctx) != T_TRUE) {
      flint_printf("Not the right displacement\n");
      gr_mat_clear(A, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return GR_TEST_FAIL;
    }
    gr_mat_clear(A, ctx);
    gr_mat_clear(D, ctx);
    gr_ctx_clear(ctx);
    i++;
  }

  flint_rand_clear(state);
  return GR_SUCCESS;
}

int test_manual_deplacement() {
  int i = 0;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  while (i < 10) {
    gr_ctx_t ctx;
    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
    flint_printf("\n:-------: Manual nxm Matrix Test :-------:\n");
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
    gr_mat_displacement(D, A, DISP_PLUS, ctx);
    flint_printf("Displacement Matrix Phi_+ (∇A):\n");
    gr_mat_print(D, ctx);
    flint_printf("\n");
    int err = gr_mat_zero(D, ctx);
    if (err) exit(FLINT_TEST_FAIL);
    gr_mat_displacement(D, A, DISP_PLUS, ctx);
    flint_printf("Displacement Matrix Phi_- (∇A):\n");
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
    gr_ctx_clear(ctx);
    i++;
  }

  flint_rand_clear(state);
  return GR_SUCCESS;
}

int test_large_matrix_deplacement_time() {
  flint_printf("\n\n:-------: Random Large nxn Matrix Execution Time Test :-------:\n");
  gr_ctx_t ctx;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
  int res = GR_SUCCESS;
  slong N = 500;
  flint_printf("Generating %d x %d matrix...\n", N, N);
  gr_mat_t A, D1, D2;
  gr_mat_init(A, N, N, ctx);
  gr_mat_init(D1, N, N, ctx);
  gr_mat_init(D2, N, N, ctx);
  clock_t start, end;
  double time_slow, time_fast;
  flint_printf("\n\nRunning Safe Method (Matrix Mul)...\n"); // matrix multiplication
  start = clock();
  gr_mat_displacement_square_safe(D1, A, ctx);
  //   flint_printf("Safe Book result D1:\n");
  //   gr_mat_print(D1, ctx);
  end = clock();
  time_slow = ((double)(end - start)) / CLOCKS_PER_SEC;
  flint_printf("\nUnoptimized Time: %f seconds\n", time_slow);
  flint_printf("Running Optimized Method (Element-wise)...\n"); // faster one
  start = clock();
  gr_mat_displacement(D2, A, DISP_PLUS, ctx);
  //   flint_printf("Safe Book result D2:\n");
  //   gr_mat_print(D2, ctx);
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
  gr_ctx_clear(ctx);
  flint_rand_clear(state);
  return res;
}

int test_toeplitz_to_G_H() {
  int i = 0;
  int error = GR_SUCCESS;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  while (i < 50) {
    slong n = n_randint(state, 100);
    slong m = n_randint(state, 100);
    gr_ctx_t ctx;
    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
    gr_mat_t A, D, G, H;
    gr_mat_init(A, n, m, ctx);
    gr_mat_init(D, n, m, ctx);
    error = gr_mat_random_toeplitz(A, state, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Error cant create Toeplitz\n");
      gr_mat_clear(A, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(G, ctx);
      gr_mat_clear(H, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_G_H(G, H, A, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Error cant create G & H\n");
      gr_mat_clear(A, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(G, ctx);
      gr_mat_clear(H, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_reconstruct_A(D, G, H, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Error cant reconstruct A\n");
      gr_mat_clear(A, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(G, ctx);
      gr_mat_clear(H, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    if (gr_mat_equal(A, D, ctx) != T_TRUE) {
      flint_printf("A & D are not equal\n");
      gr_mat_clear(A, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(G, ctx);
      gr_mat_clear(H, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return GR_TEST_FAIL;
    }
    gr_mat_clear(A, ctx);
    gr_mat_clear(D, ctx);
    gr_mat_clear(G, ctx);
    gr_mat_clear(H, ctx);
    gr_ctx_clear(ctx);
    i++;
  }
  flint_rand_clear(state);
  return error;
}

int test_quasi_toeplitz_to_G_H() {
  int i = 0;
  int error = GR_SUCCESS;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  while (i < 50) {
    slong n = n_randint(state, 100);
    slong m = n_randint(state, 100);
    gr_ctx_t ctx;
    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
    gr_mat_t A, D, G, H;
    gr_mat_init(A, n, m, ctx);
    gr_mat_init(D, n, m, ctx);
    error = gr_mat_random_quasi_toepitz(A, state, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Error cant create Toeplitz\n");
      gr_mat_clear(A, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(G, ctx);
      gr_mat_clear(H, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_G_H(G, H, A, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Error cant create G & H\n");
      gr_mat_clear(A, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(G, ctx);
      gr_mat_clear(H, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_reconstruct_A(D, G, H, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Error cant reconstruct A\n");
      gr_mat_clear(A, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(G, ctx);
      gr_mat_clear(H, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    if (gr_mat_equal(A, D, ctx) != T_TRUE) {
      flint_printf("A & D are not equal\n");
      gr_mat_clear(A, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(G, ctx);
      gr_mat_clear(H, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return GR_TEST_FAIL;
    }
    gr_mat_clear(A, ctx);
    gr_mat_clear(D, ctx);
    gr_mat_clear(G, ctx);
    gr_mat_clear(H, ctx);
    gr_ctx_clear(ctx);
    i++;
  }
  flint_rand_clear(state);
  return error;
}

int test_G_H_2x2_execution() {
  int status = GR_SUCCESS;
  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, 65537);

  // invalid case
  gr_mat_t TEMP, t1, t2;
  gr_mat_init(TEMP, 2, 2, ctx);
  status |= gr_mat_zero(TEMP, ctx);
  if (status != GR_SUCCESS) status = GR_TEST_FAIL;
  status |= gr_mat_G_H(t1, t2, TEMP, DISP_PLUS, ctx);
  if (status == GR_UNABLE) status = GR_SUCCESS;

  // valid case
  gr_mat_t A, G, H;
  gr_mat_init(A, 2, 2, ctx);
  status |= gr_mat_one(A, ctx);
  status |= gr_set_ui(gr_mat_entry_ptr(A, 0, 1, ctx), 5, ctx);
  status |= gr_mat_G_H(G, H, A, DISP_PLUS, ctx);

  gr_mat_clear(A, ctx);
  gr_mat_clear(G, ctx);
  gr_mat_clear(H, ctx);
  gr_mat_clear(TEMP, ctx);
  gr_mat_clear(t1, ctx);
  gr_mat_clear(t2, ctx);
  gr_ctx_clear(ctx);
  return status;
}

int test_toeplitz_reconstruction() {
  int i = 0;
  int res = GR_SUCCESS;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  while (i < 10) {
    flint_printf("\n\n:-------: Toeplitz Reconstruction Test :-------:\n");
    gr_ctx_t ctx;

    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

    slong n = 5;
    gr_mat_t A, A_ref, G, H, D_orig, D_from_GH, HT, B, Diff;
    gr_mat_init(A, n, n, ctx);
    gr_mat_init(A_ref, n, n, ctx);
    gr_mat_init(D_orig, n, n, ctx);
    gr_mat_init(D_from_GH, n, n, ctx);
    gr_mat_init(B, n, n, ctx);
    gr_mat_init(Diff, n, n, ctx);

    if (gr_mat_random_toeplitz(A, state, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] gr_mat_random_toeplitz failed\n");
      res = GR_TEST_FAIL;
    }
    if (gr_mat_set(A_ref, A, ctx) != GR_SUCCESS) { res = GR_TEST_FAIL; }

    flint_printf("Original Toeplitz matrix A:\n");
    gr_mat_print(A, ctx);
    flint_printf("\n");

    if (gr_mat_displacement(D_orig, A, DISP_PLUS, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] displacement failed\n");
      res = GR_TEST_FAIL;
    }
    flint_printf("True displacement D = A - Z A Z^T:\n");
    gr_mat_print(D_orig, ctx);
    flint_printf("\n");

    if (gr_mat_G_H(G, H, A, DISP_PLUS, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] gr_mat_G_H failed\n");
      res = GR_TEST_FAIL;
    }
    flint_printf("Generator G (%ld x %ld):\n", gr_mat_nrows(G, ctx), gr_mat_ncols(G, ctx));
    gr_mat_print(G, ctx);
    flint_printf("\n");
    flint_printf("Generator H (%ld x %ld):\n", gr_mat_nrows(H, ctx), gr_mat_ncols(H, ctx));
    gr_mat_print(H, ctx);
    flint_printf("\n");

    if (gr_mat_reconstruct_A(B, G, H, DISP_PLUS, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] gr_mat_reconstruct_A_safe failed\n");
      res = GR_TEST_FAIL;
    }
    flint_printf("Reconstructed matrix B:\n");
    gr_mat_print(B, ctx);
    flint_printf("\n");

    if (gr_mat_sub(Diff, B, A_ref, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] gr_mat_sub failed\n");
      res = GR_TEST_FAIL;
    }
    flint_printf("Difference B - A (should be zero):\n");
    gr_mat_print(Diff, ctx);
    flint_printf("\n");

    // Verify D = G * H^T
    gr_mat_init(HT, gr_mat_ncols(H, ctx), gr_mat_nrows(H, ctx), ctx);
    if (gr_mat_transpose(HT, H, ctx) != GR_SUCCESS || gr_mat_mul(D_from_GH, G, HT, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] transpose/mul for displacement check failed\n");
      gr_mat_clear(HT, ctx);
      res = GR_TEST_FAIL;
    }
    gr_mat_clear(HT, ctx);

    flint_printf("Displacement rebuilt from G * H^T:\n");
    gr_mat_print(D_from_GH, ctx);
    flint_printf("\n");

    if (gr_mat_sub(Diff, D_orig, D_from_GH, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] gr_mat_sub failed\n");
      res = GR_TEST_FAIL;
    }
    flint_printf("Difference (true D - D from G*H^T):\n");
    gr_mat_print(Diff, ctx);
    flint_printf("\n");

    if (gr_mat_is_zero(Diff, ctx) == T_TRUE)
      flint_printf("[GOOD] G and H correctly represent the displacement D.\n");
    else {
      flint_printf("[BAD] G*H^T does NOT equal the true displacement D!\n");
      res = GR_TEST_FAIL;
    }

    if (gr_mat_equal(B, A_ref, ctx) == T_TRUE)
      flint_printf("[SUCCESS] Reconstruction matches original A.\n");
    else {
      flint_printf("[FAILURE] Reconstruction does NOT match A.\n");
      res = GR_TEST_FAIL;
    }

    gr_mat_clear(A, ctx);
    gr_mat_clear(A_ref, ctx);
    gr_mat_clear(G, ctx);
    gr_mat_clear(H, ctx);
    gr_mat_clear(D_orig, ctx);
    gr_mat_clear(D_from_GH, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(Diff, ctx);
    gr_ctx_clear(ctx);
    i++;
  }
  flint_rand_clear(state);
  return res;
}

int test_quasi_toeplitz_reconstruction() {
  int i = 0;
  int res = GR_SUCCESS;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  while (i < 10) {
    flint_printf("\n\n:-------: Quasi-Toeplitz Full Round-Trip Reconstruction Test :-------:\n");
    gr_ctx_t ctx;

    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

    slong m = 8, n = 8;
    gr_mat_t A, A_ref, G, H, D_ref, D_from_GH, HT, B, Diff;
    gr_mat_init(A, m, n, ctx);
    gr_mat_init(A_ref, m, n, ctx);
    gr_mat_init(D_ref, m, n, ctx);
    gr_mat_init(D_from_GH, m, n, ctx);
    gr_mat_init(B, m, n, ctx);
    gr_mat_init(Diff, m, n, ctx);

    if (gr_mat_random_quasi_toepitz(A, state, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] random_quasi_toeplitz failed\n");
      res = GR_TEST_FAIL;
    }
    if (gr_mat_set(A_ref, A, ctx) != GR_SUCCESS) { res = GR_TEST_FAIL; }

    flint_printf("Original quasi-Toeplitz A:\n");
    gr_mat_print(A, ctx);
    flint_printf("\n");

    if (gr_mat_G_H(G, H, A, DISP_PLUS, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] gr_mat_G_H failed\n");
      res = GR_TEST_FAIL;
    }
    flint_printf("Generator G (%ld x %ld):\n", gr_mat_nrows(G, ctx), gr_mat_ncols(G, ctx));
    gr_mat_print(G, ctx);
    flint_printf("\nGenerator H (%ld x %ld):\n", gr_mat_nrows(H, ctx), gr_mat_ncols(H, ctx));
    gr_mat_print(H, ctx);
    flint_printf("\n");

    // Verify G * H^T == displacement(A)
    if (gr_mat_displacement(D_ref, A_ref, DISP_PLUS, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] displacement failed\n");
      res = GR_TEST_FAIL;
    }

    gr_mat_init(HT, gr_mat_ncols(H, ctx), gr_mat_nrows(H, ctx), ctx);
    if (gr_mat_transpose(HT, H, ctx) != GR_SUCCESS || gr_mat_mul(D_from_GH, G, HT, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] transpose/mul failed\n");
      gr_mat_clear(HT, ctx);
      res = GR_TEST_FAIL;
    }
    gr_mat_clear(HT, ctx);

    if (gr_mat_sub(Diff, D_ref, D_from_GH, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] gr_mat_sub failed\n");
      res = GR_TEST_FAIL;
    }
    if (gr_mat_is_zero(Diff, ctx) == T_TRUE)
      flint_printf("[OK] G * H^T == displacement(A)\n");
    else {
      flint_printf("[FAIL] G * H^T != displacement(A) — generator extraction broken!\n");
      res = GR_TEST_FAIL;
    }

    // Full round-trip: reconstruct A and compare
    if (gr_mat_reconstruct_A(B, G, H, DISP_PLUS, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] gr_mat_reconstruct_A_safe failed\n");
      res = GR_TEST_FAIL;
    }

    flint_printf("\nReconstructed B:\n");
    gr_mat_print(B, ctx);
    flint_printf("\nOriginal A (reference):\n");
    gr_mat_print(A_ref, ctx);

    if (gr_mat_sub(Diff, B, A_ref, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] gr_mat_sub failed\n");
      res = GR_TEST_FAIL;
    }
    flint_printf("\nDifference B - A (should be zero):\n");
    gr_mat_print(Diff, ctx);
    flint_printf("\n");

    if (gr_mat_equal(B, A_ref, ctx) == T_TRUE)
      flint_printf("[SUCCESS] Reconstructed B == original A.\n");
    else {
      flint_printf("[FAILURE] Reconstructed B != original A.\n");
      res = GR_TEST_FAIL;
    }

    gr_mat_clear(A, ctx);
    gr_mat_clear(A_ref, ctx);
    gr_mat_clear(G, ctx);
    gr_mat_clear(H, ctx);
    gr_mat_clear(D_ref, ctx);
    gr_mat_clear(D_from_GH, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(Diff, ctx);

    gr_ctx_clear(ctx);
    i++;
  }
  flint_rand_clear(state);
  return res;
}

int test_G_H_zero_matrix() {
  int status = GR_SUCCESS;
  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, 65537);

  slong dimensions[][2] = {{1, 1}, {2, 2}, {5, 5}, {10, 10}, {3, 5}, {5, 3}, {1, 7}, {7, 1}, {4, 9}, {9, 4}};
  slong nb_dimensions = 10;

  for (slong s = 0; s < nb_dimensions; s++) {
    if (status != GR_SUCCESS) break;

    slong rows = dimensions[s][0];
    slong cols = dimensions[s][1];

    gr_mat_t A, G, H, Reconstructed;
    gr_mat_init(A, rows, cols, ctx);
    gr_mat_init(Reconstructed, rows, cols, ctx);
    status |= gr_mat_zero(A, ctx);

    int gh_status = gr_mat_G_H(G, H, A, DISP_PLUS, ctx);
    if (gh_status != GR_SUCCESS) {
      status = GR_TEST_FAIL;
      gr_mat_clear(A, ctx);
      gr_mat_clear(Reconstructed, ctx);
      continue;
    }

    if (gr_mat_ncols(G, ctx) < 1 || gr_mat_ncols(H, ctx) < 1) {
      status = GR_TEST_FAIL;
      gr_mat_clear(A, ctx);
      gr_mat_clear(G, ctx);
      gr_mat_clear(H, ctx);
      gr_mat_clear(Reconstructed, ctx);
      continue;
    }

    // G should be rows x rank and H should be cols x rank
    if (gr_mat_nrows(G, ctx) != rows) status = GR_TEST_FAIL;
    if (gr_mat_nrows(H, ctx) != cols) status = GR_TEST_FAIL;

    if (gr_mat_is_zero(G, ctx) == T_FALSE) status = GR_TEST_FAIL;
    if (gr_mat_is_zero(H, ctx) == T_FALSE) status = GR_TEST_FAIL;

    // big verification
    int rec_status = gr_mat_reconstruct_A(Reconstructed, G, H, DISP_PLUS, ctx);
    if (rec_status != GR_SUCCESS)
      status = GR_TEST_FAIL;
    else if (gr_mat_equal(Reconstructed, A, ctx) == T_FALSE) // A = disp(A)
      status = GR_TEST_FAIL;

    gr_mat_clear(A, ctx);
    gr_mat_clear(G, ctx);
    gr_mat_clear(H, ctx);
    gr_mat_clear(Reconstructed, ctx);
  }
  gr_ctx_clear(ctx);
  return status;
}



int test_toeplitz_reconstruction_v2() {
  int i = 0;
  int res = GR_SUCCESS;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  while (i < 10) {
     gr_ctx_t ctx;

    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

    slong n = 5;
    gr_mat_t A, G, H, B;
    gr_mat_init(A, n, n, ctx);
    gr_mat_init(B, n, n, ctx);

    if (gr_mat_random_toeplitz(A, state, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] gr_mat_random_toeplitz failed\n");
      res = GR_TEST_FAIL;
    }

    if (gr_mat_G_H(G, H, A, DISP_PLUS, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] gr_mat_G_H failed\n");
      res = GR_TEST_FAIL;
    }

    if (gr_mat_reconstruct_A_v2(B, G, H, DISP_PLUS, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] gr_mat_reconstruct_A_safe failed\n");
      res = GR_TEST_FAIL;
    }
    if (gr_mat_equal(B, A, ctx) == T_TRUE);
    else {
      flint_printf("[FAILURE] Reconstructed B != original A.\n");
      flint_printf("Original A:\n");
      gr_mat_print(A, ctx);
      flint_printf("\n");
      flint_printf("Original G:\n");
      gr_mat_print(G, ctx);
      flint_printf("\n");
      flint_printf("Original H:\n");
      gr_mat_print(H, ctx);
      flint_printf("\n");
      flint_printf("Reconstructed B:\n");
      gr_mat_print(B, ctx);
      flint_printf("\n");
      res = GR_TEST_FAIL;
      i = 9999;
    }

    gr_mat_clear(A, ctx);
    gr_mat_clear(G, ctx);
    gr_mat_clear(H, ctx);
    gr_mat_clear(B, ctx);
    gr_ctx_clear(ctx);
    i++;
  }
  flint_rand_clear(state);
  return res;
}

int test_quasi_toeplitz_reconstruction_v2() {
  int i = 0;
  int res = GR_SUCCESS;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  while (i < 10) {
    flint_printf("\n\n:-------: Quasi-Toeplitz Full Round-Trip Reconstruction Test :-------:\n");
    gr_ctx_t ctx;

    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

    slong m = 8, n = 8;
    gr_mat_t A, G, H, B;
    gr_mat_init(A, m, n, ctx);


    gr_mat_init(B, m, n, ctx);


    if (gr_mat_random_quasi_toepitz(A, state, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] random_quasi_toeplitz failed\n");
      res = GR_TEST_FAIL;
    }
    if (gr_mat_G_H(G, H, A, DISP_PLUS, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] gr_mat_G_H failed\n");
      res = GR_TEST_FAIL;
    }
    if (gr_mat_reconstruct_A_v2(B, G, H, DISP_PLUS, ctx) != GR_SUCCESS) {
      flint_printf("[ERROR] gr_mat_reconstruct_A_safe failed\n");
      res = GR_TEST_FAIL;
    }

    if (gr_mat_equal(B, A, ctx) == T_TRUE); 
    else {
      flint_printf("[FAILURE] Reconstructed B != original A.\n");
      flint_printf("Original A:\n");
      gr_mat_print(A, ctx);
      flint_printf("\n");
      flint_printf("Original G:\n");
      gr_mat_print(G, ctx);
      flint_printf("\n");
      flint_printf("Original H:\n");
      gr_mat_print(H, ctx);
      flint_printf("\n");
      flint_printf("Reconstructed B:\n");
      gr_mat_print(B, ctx);
      flint_printf("\n");
      res = GR_TEST_FAIL;
      i = 9999;
    }
    gr_mat_clear(A, ctx);
    gr_mat_clear(G, ctx);
    gr_mat_clear(H, ctx);
    gr_mat_clear(B, ctx);
    gr_ctx_clear(ctx);
    i++;
  }
  flint_rand_clear(state);
  return res;
}

void usage(char *argv[]) {
  fprintf(stderr, "Usage: %s <test_name>\n", argv[0]);
  fprintf(stderr, "Available tests:\n");
  fprintf(stderr, "  - toeplitz_deplacement\n");
  fprintf(stderr, "  - manual_deplacement\n");
  fprintf(stderr, "  - large_deplacement\n");
  fprintf(stderr, "  - toeplitz_G_H\n");
  fprintf(stderr, "  - toeplitz_reconstruction\n");
  fprintf(stderr, "  - quasi_toeplitz_reconstruction\n");
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    usage(argv);
    return GR_UNABLE;
  }

  // start test
  fprintf(stderr, "=> Start test \"%s\"\n", argv[1]);
  int ok = GR_SUCCESS;
  if (strcmp("toeplitz_deplacement", argv[1]) == 0) {
    ok = test_toeplitz_deplacement();
  } else if (strcmp("manual_deplacement", argv[1]) == 0) {
    ok = test_manual_deplacement();
  } else if (strcmp("large_deplacement", argv[1]) == 0) {
    ok = test_large_matrix_deplacement_time();
  } else if (strcmp("toeplitz_G_H", argv[1]) == 0) {
    ok = test_toeplitz_to_G_H();
  } else if (strcmp("quasi_toeplitz_G_H", argv[1]) == 0) {
    ok = test_quasi_toeplitz_to_G_H();
  } else if (strcmp("toeplitz_reconstruction", argv[1]) == 0) {
    ok = test_toeplitz_reconstruction();
  } else if (strcmp("quasi_toeplitz_reconstruction", argv[1]) == 0) {
    ok = test_quasi_toeplitz_reconstruction();
  } else if (strcmp("G_H_2x2_execution", argv[1]) == 0) {
    ok = test_G_H_2x2_execution();
  } else if (strcmp("G_H_zero_matrix", argv[1]) == 0) {
    ok = test_G_H_zero_matrix();
  } else if (strcmp("toeplitz_reconstruction_v2", argv[1]) == 0) {
    ok = test_toeplitz_reconstruction_v2();
  } else if (strcmp("quasi_toeplitz_reconstruction_v2", argv[1]) == 0) {
    ok = test_quasi_toeplitz_reconstruction_v2();
  } else {
    fprintf(stderr, "Error: test \"%s\" not found!\n", argv[1]);
    exit(EXIT_FAILURE);
  }

  if (!ok) {
    fprintf(stderr, "Test \"%s\" finished: SUCCESS\n", argv[1]);
    return EXIT_SUCCESS;
  } else {
    fprintf(stderr, "Test \"%s\" finished: FAILURE\n", argv[1]);
    return EXIT_FAILURE;
  }
}
