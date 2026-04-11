#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"
#include "flint/ulong_extras.h"
#include "random_toeplitz.h"
#include "inverse_toeplitz.h"

int test_inverse_base_case() {
  int i = 0;
  int status = GR_SUCCESS;
  flint_rand_t state;
  flint_rand_init(state);

  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

  slong test_rank = 3;

  gr_mat_t A, A_inv, G_A, H_A, G_D, H_D, Check;
  gr_mat_init(A, 1, 1, ctx);
  gr_mat_init(A_inv, 1, 1, ctx);
  
  gr_mat_init(G_A, 1, test_rank, ctx);
  gr_mat_init(H_A, 1, test_rank, ctx);
  
  gr_mat_init(G_D, 1, 1, ctx);
  gr_mat_init(H_D, 1, 1, ctx);
  gr_mat_init(Check, 1, 1, ctx);

  while (i < 50) {
    status |= gr_mat_randtest(G_A, state, ctx);
    status |= gr_mat_randtest(H_A, state, ctx);

    status |= gr_mat_reconstruct_A(A, G_A, H_A, DISP_PLUS, ctx);

    if (gr_mat_is_zero(A, ctx) == T_TRUE) continue;

    status |= gr_mat_inverse_toeplitz(G_D, H_D, G_A, H_A, ctx);
    status |= gr_mat_reconstruct_A(A_inv, G_D, H_D, DISP_PLUS, ctx);
    status |= gr_mat_mul(Check, A, A_inv, ctx);

    if (gr_mat_is_one(Check, ctx) != T_TRUE) {
      status = GR_TEST_FAIL;
      break;
    }
    i++;
  }

  gr_mat_clear(A, ctx);
  gr_mat_clear(A_inv, ctx);
  gr_mat_clear(G_A, ctx);
  gr_mat_clear(H_A, ctx);
  gr_mat_clear(G_D, ctx);
  gr_mat_clear(H_D, ctx);
  gr_mat_clear(Check, ctx);
  gr_ctx_clear(ctx);
  flint_rand_clear(state);
  return status;
}

int test_inverse_2x2() {
  int status = GR_SUCCESS;
  flint_rand_t state;
  flint_rand_init(state);
  gr_ctx_t ctx;
  // gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
  gr_ctx_init_nmod(ctx, 11);

  gr_mat_t A, A_inv, Check;
  gr_mat_t G_A, H_A, G_D, H_D;

  gr_mat_init(A, 2, 2, ctx);
  status |= gr_set_ui(gr_mat_entry_ptr(A, 0, 0, ctx), 1, ctx);
  status |= gr_set_ui(gr_mat_entry_ptr(A, 0, 1, ctx), 3, ctx);
  status |= gr_set_ui(gr_mat_entry_ptr(A, 1, 0, ctx), 2, ctx);
  status |= gr_set_ui(gr_mat_entry_ptr(A, 1, 1, ctx), 4, ctx);

  printf("A =\n");
  gr_mat_print(A, ctx);
  printf("\n");

  gr_mat_init(G_A, 2, 4, ctx);
  gr_mat_init(H_A, 2, 4, ctx);

  status = gr_mat_G_H(G_A, H_A, A, DISP_PLUS, ctx);
  if (status != GR_SUCCESS) goto cleanup;

  printf("G_A:\n");
  gr_mat_print(G_A, ctx);
  printf("\nH_A:\n");
  gr_mat_print(H_A, ctx);
  printf("\n");

  gr_mat_init(G_D, 1, 1, ctx);
  gr_mat_init(H_D, 1, 1, ctx);

  status = gr_mat_inverse_toeplitz(G_D, H_D, G_A, H_A, ctx);

  if (status != GR_SUCCESS) {
    printf("gr_mat_inverse_toeplitz failed with status = %d\n", status);
    goto cleanup;
  }

  printf("G_D (generators of A^-1):\n");
  gr_mat_print(G_D, ctx);
  printf("\nH_D (generators of A^-1):\n");
  gr_mat_print(H_D, ctx);
  printf("\n");

  // Reconstruct inverse
  gr_mat_init(A_inv, 2, 2, ctx);
  status = gr_mat_reconstruct_A(A_inv, G_D, H_D, DISP_PLUS, ctx);
  if (status != GR_SUCCESS) goto cleanup;

  printf("A_inv (reconstructed):\n");
  gr_mat_print(A_inv, ctx);
  printf("\n");

  gr_mat_init(Check, 2, 2, ctx);
  status = gr_mat_mul(Check, A, A_inv, ctx);
  if (status != GR_SUCCESS) goto cleanup;

  printf("Check = A * A_inv =\n");
  gr_mat_print(Check, ctx);
  printf("\n");

  if (gr_mat_is_one(Check, ctx) != T_TRUE) {
    printf("Product is not identity\n");
    status = GR_TEST_FAIL;
  }

cleanup:
  gr_mat_clear(A, ctx);
  gr_mat_clear(A_inv, ctx);
  gr_mat_clear(Check, ctx);
  gr_mat_clear(G_A, ctx);
  gr_mat_clear(H_A, ctx);
  gr_mat_clear(G_D, ctx);
  gr_mat_clear(H_D, ctx);
  gr_ctx_clear(ctx);
  flint_rand_clear(state);
  return status;
}

int test_inverse_3x3() {
  int status = GR_SUCCESS;
  flint_rand_t state;
  flint_rand_init(state);
  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, 11);

  gr_mat_t A, A_inv, A_inv_ref, Check;
  gr_mat_t G_A, H_A, G_D, H_D;

  gr_mat_init(A, 3, 3, ctx);
  gr_mat_init(A_inv, 3, 3, ctx);
  gr_mat_init(A_inv_ref, 3, 3, ctx);
  gr_mat_init(Check, 3, 3, ctx);
  gr_mat_init(G_A, 3, 6, ctx);
  gr_mat_init(H_A, 3, 6, ctx);
  gr_mat_init(G_D, 1, 1, ctx);
  gr_mat_init(H_D, 1, 1, ctx);

  status |= gr_set_ui(gr_mat_entry_ptr(A, 0, 0, ctx), 1, ctx);
  status |= gr_set_ui(gr_mat_entry_ptr(A, 0, 1, ctx), 2, ctx);
  status |= gr_set_ui(gr_mat_entry_ptr(A, 0, 2, ctx), 3, ctx);
  status |= gr_set_ui(gr_mat_entry_ptr(A, 1, 0, ctx), 4, ctx);
  status |= gr_set_ui(gr_mat_entry_ptr(A, 1, 1, ctx), 5, ctx);
  status |= gr_set_ui(gr_mat_entry_ptr(A, 1, 2, ctx), 6, ctx);
  status |= gr_set_ui(gr_mat_entry_ptr(A, 2, 0, ctx), 7, ctx);
  status |= gr_set_ui(gr_mat_entry_ptr(A, 2, 1, ctx), 8, ctx);
  status |= gr_set_ui(gr_mat_entry_ptr(A, 2, 2, ctx), 10, ctx);

  printf("A =\n");
  gr_mat_print(A, ctx);
  printf("\n");

  if (status == GR_SUCCESS) {
    status = gr_mat_G_H(G_A, H_A, A, DISP_PLUS, ctx);
    printf("G_A:\n");
    gr_mat_print(G_A, ctx);
    printf("\nH_A:\n");
    gr_mat_print(H_A, ctx);
    printf("\n");
  }

  if (status == GR_SUCCESS) {
    status = gr_mat_inverse_toeplitz(G_D, H_D, G_A, H_A, ctx);
    if (status != GR_SUCCESS) printf("gr_mat_inverse_toeplitz failed with status = %d\n", status);
  }

  if (status == GR_SUCCESS) {
    printf("G_D (generators of A^-1):\n");
    gr_mat_print(G_D, ctx);
    printf("\nH_D (generators of A^-1):\n");
    gr_mat_print(H_D, ctx);
    printf("\n");

    status = gr_mat_reconstruct_A(A_inv, G_D, H_D, DISP_PLUS, ctx);
    if (status != GR_SUCCESS) printf("reconstruct failed\n");
  }

  if (status == GR_SUCCESS) {
    printf("A_inv (reconstructed):\n");
    gr_mat_print(A_inv, ctx);
    printf("\n");

    status = gr_mat_inv(A_inv_ref, A, ctx);
    if (status != GR_SUCCESS) printf("gr_mat_inv reference failed (singular?)\n");
  }

  if (status == GR_SUCCESS) {
    printf("A_inv reference (gr_mat_inv):\n");
    gr_mat_print(A_inv_ref, ctx);
    printf("\n");

    status = gr_mat_mul(Check, A, A_inv, ctx);
  }

  if (status == GR_SUCCESS) {
    printf("Check = A * A_inv =\n");
    gr_mat_print(Check, ctx);
    printf("\n");

    if (gr_mat_is_one(Check, ctx) != T_TRUE) {
      printf("Product is not identity\n");
      status = GR_TEST_FAIL;
    } else {
      printf("SUCCESS: A * A_inv = I\n");
    }
  }

  gr_mat_clear(A, ctx);
  gr_mat_clear(A_inv, ctx);
  gr_mat_clear(A_inv_ref, ctx);
  gr_mat_clear(Check, ctx);
  gr_mat_clear(G_A, ctx);
  gr_mat_clear(H_A, ctx);
  gr_mat_clear(G_D, ctx);
  gr_mat_clear(H_D, ctx);
  gr_ctx_clear(ctx);
  flint_rand_clear(state);
  return status;
}

int test_inverse_full() {
  int i = 0;
  int status = GR_SUCCESS;
  flint_rand_t state;
  flint_rand_init(state);
  int size = 5;
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);

  while (i < 10) {
    gr_mat_t A, B, G, H, T, U, C;
    gr_ctx_t ctx;
    gr_ptr det;

    slong n = n_randint(state, size) + 1;
    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

    gr_mat_init(A, n, n, ctx);
    gr_mat_init(B, n, n, ctx);
    gr_mat_init(C, n, n, ctx);
    GR_TMP_INIT(det, ctx);

    do {
      status |= gr_mat_random_toeplitz(A, state, ctx);
      if (status != GR_SUCCESS) goto fail;
      status |= gr_mat_det(det, A, ctx);
    } while (gr_is_zero(det, ctx) == T_TRUE);

    gr_mat_init(G, n, 2, ctx);
    gr_mat_init(H, n, 2, ctx);
    gr_mat_init(T, n, 2, ctx);
    gr_mat_init(U, n, 2, ctx);

    status |= gr_mat_G_H(G, H, A, DISP_PLUS, ctx);
    if (status == GR_SUCCESS) status |= gr_mat_inverse_toeplitz(T, U, G, H, ctx);

    if (status == GR_SUCCESS) status |= gr_mat_reconstruct_A(B, T, U, DISP_PLUS, ctx);

    if (status == GR_SUCCESS) {
      status |= gr_mat_mul(C, A, B, ctx);
      if (gr_mat_is_one(C, ctx) != T_TRUE) {
        printf("Failure: A * A_inv != I at size %ld\n", n);
        status = GR_TEST_FAIL;
      }
    }

  fail:
    GR_TMP_CLEAR(det, ctx);
    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);
    gr_mat_clear(G, ctx);
    gr_mat_clear(H, ctx);
    gr_mat_clear(T, ctx);
    gr_mat_clear(U, ctx);
    gr_ctx_clear(ctx);

    if (status != GR_SUCCESS) break;

    i++;
    size *= 2;
  }
  flint_rand_clear(state);
  return status;
}

void usage(char *argv[]) { fprintf(stderr, "Usage: %s <test_name>\n", argv[0]); }

int main(int argc, char *argv[]) {
  if (argc < 2) {
    usage(argv);
    return GR_UNABLE;
  }
  // start test
  fprintf(stderr, "=> Start test \"%s\"\n", argv[1]);
  int ok = GR_SUCCESS;
  if (strcmp("inverse_base_case", argv[1]) == 0) {
    ok = test_inverse_base_case();
  } else if (strcmp("inverse_2x2", argv[1]) == 0) {
    ok = test_inverse_2x2();
  } else if (strcmp("inverse_3x3", argv[1]) == 0) {
    ok = test_inverse_3x3();
  } else if (strcmp("inverse_full", argv[1]) == 0) {
    ok = test_inverse_full();
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
