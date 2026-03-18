#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"
#include "flint/ulong_extras.h"
#include "toeplitz_inverse.h"

int test_inverse_base_case() {
  int i = 0;
  int status = GR_SUCCESS;
  flint_rand_t state;
  flint_rand_init(state);

  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

  gr_mat_t A, A_inv, G_A, H_A, G_D, H_D, Check;
  gr_mat_init(A, 1, 1, ctx);
  gr_mat_init(A_inv, 1, 1, ctx);
  gr_mat_init(G_A, 1, 1, ctx);
  gr_mat_init(H_A, 1, 1, ctx);
  gr_mat_init(G_D, 1, 1, ctx);
  gr_mat_init(H_D, 1, 1, ctx);
  gr_mat_init(Check, 1, 1, ctx);

  while (i < 50) {
    status |= gr_mat_zero(A, ctx);
    do { status |= gr_mat_randtest(A, state, ctx); } while (gr_mat_is_zero(A, ctx) == T_TRUE);
    status |= gr_mat_G_H(G_A, H_A, A, DISP_PLUS, ctx);
    status |= gr_toeplitz_inverse(G_D, H_D, G_A, H_A, ctx);
    status |= gr_mat_reconstruct_A(A_inv, G_D, H_D, DISP_PLUS, ctx);
    status |= gr_mat_mul(Check, A, A_inv, ctx);

    if (gr_mat_is_one(Check, ctx) != T_TRUE) {
      printf("Failed at the iteration %d\n", i);
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
  int i = 0;
  int status = GR_SUCCESS;
  flint_rand_t state;
  flint_rand_init(state);

  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

  gr_mat_t A, A_inv, G_A, H_A, G_D, H_D, Check;
  gr_mat_init(A, 2, 2, ctx);
  gr_mat_init(A_inv, 2, 2, ctx);
  gr_mat_init(G_A, 2, 2, ctx);
  gr_mat_init(H_A, 2, 2, ctx);
  gr_mat_init(G_D, 2, 2, ctx);
  gr_mat_init(H_D, 2, 2, ctx);
  gr_mat_init(Check, 2, 2, ctx);

  gr_ptr d = gr_heap_init(ctx);

  while (i < 20) {
    do {
      status |= gr_mat_randtest(A, state, ctx);
      status |= gr_mat_det(d, A, ctx);
    } while (gr_is_zero(d, ctx) == T_TRUE);
    status |= gr_mat_G_H(G_A, H_A, A, DISP_PLUS, ctx);
    status |= gr_toeplitz_inverse(G_D, H_D, G_A, H_A, ctx);
    status |= gr_mat_reconstruct_A(A_inv, G_D, H_D, DISP_PLUS, ctx);
    status |= gr_mat_mul(Check, A, A_inv, ctx);

    if (gr_mat_is_one(Check, ctx) != T_TRUE) {
      printf("Failed at the iteration %d\n", i);
      status = GR_TEST_FAIL;
      break;
    }
    i++;
  }

  // Cleanup
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
