#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/ulong_extras.h"
#include "matrix_aux.h"
#include "multiplication.h"
#include "random_toeplitz.h"

int test_multiplication_toeplitz() {
  int i = 0;
  int error = GR_SUCCESS;
  flint_rand_t state;
  gr_ctx_t ctx;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
  while (i < 50) {
    gr_mat_t A, B, C, D, G_a, H_a, G_b, H_b, G_c, H_c;
    slong n = n_randint(state, 15);
    slong m = n_randint(state, 15);
    slong k = n_randint(state, 15);

    gr_mat_init(A, n, m, ctx);
    gr_mat_init(B, m, k, ctx);
    error = gr_mat_random_toeplitz(A, state, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to create a Toeplitz A, FIX IT FAST\n");
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_random_toeplitz(B, state, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to create a Toeplitz A, FIX IT FAST\n");
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }

    gr_mat_init(C, gr_mat_nrows(A, ctx), gr_mat_ncols(B, ctx), ctx);
    gr_mat_init(D, gr_mat_nrows(A, ctx), gr_mat_ncols(B, ctx), ctx);
    error = gr_mat_G_H(G_a, H_a, A, DISP_PLUS, ctx);

    if (error != GR_SUCCESS) {
      flint_printf("Failed to generate G and H of A, FIX IT FAST\n");
      gr_mat_clear(C, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }

    error = gr_mat_G_H(G_b, H_b, B, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to generate G and H of B, FIX IT FAST\n");
      gr_mat_clear(C, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(G_b, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(H_b, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_mul(C, A, B, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to multiply with flint, size ??\n");
      gr_mat_clear(C, ctx);
      gr_mat_clear(A, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(G_b, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(H_b, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_mul_generator(G_c, H_c, G_a, H_a, G_b, H_b, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to multiply, FIX IT FAST\n");
      gr_mat_clear(C, ctx);
      gr_mat_clear(A, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(G_b, ctx);
      gr_mat_clear(G_c, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(H_b, ctx);
      gr_mat_clear(H_c, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_reconstruct_A(D, G_c, H_c, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to reconstruct D, FIX IT FAST\n");
      gr_mat_clear(C, ctx);
      gr_mat_clear(A, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(G_b, ctx);
      gr_mat_clear(G_c, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(H_b, ctx);
      gr_mat_clear(H_c, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    if (gr_mat_equal(D, C, ctx) != T_TRUE) {
      flint_printf("C != D\n");
      flint_printf("Matrice A = \n");
      gr_mat_print(A, ctx);
      flint_printf("\n");
      flint_printf("Matrice B = \n");
      gr_mat_print(B, ctx);
      flint_printf("\n");
      flint_printf("Matrice C = \n");
      gr_mat_print(C, ctx);
      flint_printf("\n");
      flint_printf("Matrice G_a = \n");
      gr_mat_print(G_a, ctx);
      flint_printf("\n");
      flint_printf("Matrice H_a = \n");
      gr_mat_print(H_a, ctx);
      flint_printf("\n");
      flint_printf("Matrice G_b = \n");
      gr_mat_print(G_b, ctx);
      flint_printf("\n");
      flint_printf("Matrice H_b = \n");
      gr_mat_print(H_b, ctx);
      flint_printf("\n");
      flint_printf("Matrice G_c = \n");
      gr_mat_print(G_c, ctx);
      flint_printf("\n");
      flint_printf("Matrice H_c = \n");
      gr_mat_print(H_c, ctx);
      flint_printf("\n");
      flint_printf("Matrice D = \n");
      gr_mat_print(D, ctx);
      flint_printf("\n");
      error = 1;
      i = 9998;
    }
    gr_mat_clear(D, ctx);
    gr_mat_clear(C, ctx);
    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(G_a, ctx);
    gr_mat_clear(G_b, ctx);
    gr_mat_clear(G_c, ctx);
    gr_mat_clear(H_a, ctx);
    gr_mat_clear(H_b, ctx);
    gr_mat_clear(H_c, ctx);
    i++;
  }
  flint_rand_clear(state);
  gr_ctx_clear(ctx);
  return error;
}

int test_multiplication_quasi_toeplitz() {
  int i = 0;
  int error = GR_SUCCESS;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  while (i < 100) {
    gr_mat_t A, B, C, D, G_a, H_a, G_b, H_b, G_c, H_c;
    gr_ctx_t ctx;
    slong n = n_randint(state, 100);
    slong m = n_randint(state, 100);
    slong k = n_randint(state, 100);
    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
    gr_mat_init(A, n, m, ctx);
    gr_mat_init(B, m, k, ctx);
    error = gr_mat_random_quasi_toepitz(A, state, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to create a quasi Toeplitz A, FIX IT FAST\n");
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_random_quasi_toepitz(B, state, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to create a quasi Toeplitz B, FIX IT FAST\n");
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    gr_mat_init(C, gr_mat_nrows(A, ctx), gr_mat_ncols(B, ctx), ctx);
    gr_mat_init(D, gr_mat_nrows(A, ctx), gr_mat_ncols(B, ctx), ctx);
    error = gr_mat_G_H(G_a, H_a, A, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to generate G and H of A, FIX IT FAST\n");
      gr_mat_clear(C, ctx);
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_G_H(G_b, H_b, B, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to generate G and H of B, FIX IT FAST\n");
      gr_mat_clear(C, ctx);
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(G_b, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(H_b, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_mul(C, A, B, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to multiply with flint, size ??\n");
      gr_mat_clear(C, ctx);
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(G_b, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(H_b, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_mul_generator(G_c, H_c, G_a, H_a, G_b, H_b, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to multiply, FIX IT FAST\n");
      gr_mat_clear(C, ctx);
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(G_b, ctx);
      gr_mat_clear(G_c, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(H_b, ctx);
      gr_mat_clear(H_c, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_reconstruct_A(D, G_c, H_c, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to reconstruct D, FIX IT FAST\n");
      gr_mat_clear(C, ctx);
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(G_b, ctx);
      gr_mat_clear(G_c, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(H_b, ctx);
      gr_mat_clear(H_c, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    if (gr_mat_equal(D, C, ctx) != T_TRUE) {
      flint_printf("C != D\n");
      flint_printf("Matrice A = \n");
      gr_mat_print(A, ctx);
      flint_printf("\n");
      flint_printf("Matrice B = \n");
      gr_mat_print(B, ctx);
      flint_printf("\n");
      flint_printf("Matrice C = \n");
      gr_mat_print(C, ctx);
      flint_printf("\n");
      flint_printf("Matrice G_a = \n");
      gr_mat_print(G_a, ctx);
      flint_printf("\n");
      flint_printf("Matrice H_a = \n");
      gr_mat_print(H_a, ctx);
      flint_printf("\n");
      flint_printf("Matrice G_b = \n");
      gr_mat_print(G_b, ctx);
      flint_printf("\n");
      flint_printf("Matrice H_b = \n");
      gr_mat_print(H_b, ctx);
      flint_printf("\n");
      flint_printf("Matrice G_c = \n");
      gr_mat_print(G_c, ctx);
      flint_printf("\n");
      flint_printf("Matrice H_c = \n");
      gr_mat_print(H_c, ctx);
      flint_printf("\n");
      flint_printf("Matrice D = \n");
      gr_mat_print(D, ctx);
      flint_printf("\n");
      error = GR_TEST_FAIL;
      i = 9998;
    }

    gr_mat_clear(C, ctx);
    gr_mat_clear(D, ctx);
    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(G_a, ctx);
    gr_mat_clear(G_b, ctx);
    gr_mat_clear(G_c, ctx);
    gr_mat_clear(H_a, ctx);
    gr_mat_clear(H_b, ctx);
    gr_mat_clear(H_c, ctx);
    gr_ctx_clear(ctx);
    i++;
  }
  flint_rand_clear(state);
  return error;
}

int test_multiplication_vector() {
  int i = 0;
  int error = GR_SUCCESS;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  while (i < 10) {
    gr_mat_t A, X, C, D, G_a, H_a;
    gr_ctx_t ctx;
    slong n = n_randint(state, 500);
    slong m = n_randint(state, 600);
    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
    gr_mat_init(A, n, m, ctx);
    gr_mat_init(X, m, 1, ctx);
    gr_mat_init(C, n, 1, ctx);
    gr_mat_init(D, n, 1, ctx);
    error = gr_mat_random_toeplitz(A, state, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to create to create a Toeplitz, FIX IT FAST\n");
      gr_mat_clear(C, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(A, ctx);
      gr_mat_clear(X, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }

    for (int j = 0; j < gr_mat_ncols(A, ctx); j++) {
      error = gr_randtest_not_zero(gr_mat_entry_ptr(X, j, 0, ctx), state, ctx);
      if (error != GR_SUCCESS) {
        flint_printf("Failed to create X, make sure you put the right size\n");
        gr_mat_clear(C, ctx);
        gr_mat_clear(D, ctx);
        gr_mat_clear(A, ctx);
        gr_mat_clear(X, ctx);
        gr_ctx_clear(ctx);
        flint_rand_clear(state);
        return error;
      }
    }

    error = gr_mat_mul(C, A, X, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to do the multiplication flint, make sur you give the right size\n");
      gr_mat_clear(C, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(A, ctx);
      gr_mat_clear(X, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }

    error = gr_mat_G_H(G_a, H_a, A, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to get the generators of A, FIX IT FAST !\n");
      gr_mat_clear(C, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(A, ctx);
      gr_mat_clear(X, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }

    error = gr_mat_mul_vector(D, G_a, H_a, X, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to do the multiplication, FIX IT FAST !\n");
      gr_mat_clear(C, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(A, ctx);
      gr_mat_clear(X, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    if (gr_mat_equal(D, C, ctx) != T_TRUE) {
      flint_printf("C != D\n");
      flint_printf("Matrice A = \n");
      gr_mat_print(A, ctx);
      flint_printf("\n");
      flint_printf("Matrice X = \n");
      gr_mat_print(X, ctx);
      flint_printf("\n");
      flint_printf("Matrice C = \n");
      gr_mat_print(C, ctx);
      flint_printf("\n");
      flint_printf("Matrice G_a = \n");
      gr_mat_print(G_a, ctx);
      flint_printf("\n");
      flint_printf("Matrice H_a = \n");
      gr_mat_print(H_a, ctx);
      flint_printf("\n");
      flint_printf("Matrice D = \n");
      gr_mat_print(D, ctx);
      flint_printf("\n");
      error = 1;
      i = 9998;
    }
    gr_mat_clear(C, ctx);
    gr_mat_clear(D, ctx);
    gr_mat_clear(A, ctx);
    gr_mat_clear(X, ctx);
    gr_mat_clear(G_a, ctx);
    gr_mat_clear(H_a, ctx);
    gr_ctx_clear(ctx);
    i++;
  }
  flint_rand_clear(state);
  return error;
}

int test_apply_Z() {
  int status = GR_SUCCESS;
  gr_ctx_t ctx;
  flint_rand_t state;

  gr_ctx_init_nmod(ctx, 11);
  flint_rand_init(state);

  slong n = 5;
  gr_mat_t V, ZV, ZtV;
  gr_mat_init(V, n, 1, ctx);
  gr_mat_init(ZV, n, 1, ctx);
  gr_mat_init(ZtV, n, 1, ctx);

  status |= gr_mat_randtest(V, state, ctx);

  status |= gr_mat_apply_Z(ZV, V, ctx);
  if (gr_is_zero(gr_mat_entry_ptr(ZV, 0, 0, ctx), ctx) != T_TRUE) { status = GR_TEST_FAIL; }
  for (slong i = 1; i < n; i++) {
    if (gr_equal(gr_mat_entry_ptr(ZV, i, 0, ctx), gr_mat_entry_ptr(V, i - 1, 0, ctx), ctx) != T_TRUE) {
      status = GR_TEST_FAIL;
    }
  }

  status |= gr_mat_apply_Zt(ZtV, V, ctx);

  // last element should be 0
  if (gr_is_zero(gr_mat_entry_ptr(ZtV, n - 1, 0, ctx), ctx) != T_TRUE) { status = GR_TEST_FAIL; }
  // ZtV[i] should be V[i+1]
  for (slong i = 0; i < n - 1; i++) {
    if (gr_equal(gr_mat_entry_ptr(ZtV, i, 0, ctx), gr_mat_entry_ptr(V, i + 1, 0, ctx), ctx) != T_TRUE) {
      printf("Failure: Zt shift-up mismatch at index %ld\n", i);
      status = GR_TEST_FAIL;
    }
  }

  gr_mat_t V_inplace;
  gr_mat_init(V_inplace, n, 1, ctx);
  status |= gr_mat_set(V_inplace, V, ctx);

  status |= gr_mat_apply_Zt(V_inplace, V_inplace, ctx);
  if (gr_mat_equal(V_inplace, ZtV, ctx) != T_TRUE) status = GR_TEST_FAIL;

  gr_mat_clear(V, ctx);
  gr_mat_clear(ZV, ctx);
  gr_mat_clear(ZtV, ctx);
  gr_mat_clear(V_inplace, ctx);
  flint_rand_clear(state);
  gr_ctx_clear(ctx);
  return status;
}

void usage(char *argv[]) {
  fprintf(stderr, "Usage: %s <test_name>\n", argv[0]);
  fprintf(stderr, "Available tests:\n");
  fprintf(stderr, "  - multiplication_toeplitz\n");
  fprintf(stderr, "  - multiplication_quasi_toeplitz\n");
  fprintf(stderr, "  - multiplication_vector\n");
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    usage(argv);
    return GR_UNABLE;
  }
  // start test
  fprintf(stderr, "=> Start test \"%s\"\n", argv[1]);
  int ok = GR_SUCCESS;
  if (strcmp("multiplication_toeplitz", argv[1]) == 0) {
    ok = test_multiplication_toeplitz();
  } else if (strcmp("multiplication_quasi_toeplitz", argv[1]) == 0) {
    ok = test_multiplication_quasi_toeplitz();
  } else if (strcmp("multiplication_vector", argv[1]) == 0) {
    ok = test_multiplication_vector();
  } else if (strcmp("apply_Z", argv[1]) == 0) {
    ok = test_apply_Z();
  } else {
    fprintf(stderr, "Error: test \"%s\" not found!\n", argv[1]);
    exit(EXIT_FAILURE);
  }
  flint_cleanup();
  if (!ok) {
    fprintf(stderr, "Test \"%s\" finished: SUCCESS\n", argv[1]);
    return EXIT_SUCCESS;
  } else {
    fprintf(stderr, "Test \"%s\" finished: FAILURE\n", argv[1]);
    return EXIT_FAILURE;
  }
}
