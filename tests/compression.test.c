#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "displacement_matrices.h"
#include "multiplication.h"
#include "compression.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"
#include "flint/ulong_extras.h"
#include "random_toeplitz.h"

int test_compression_then_multiply() {
  int i = 0;
  int n = 20;
  int m = n;
  int status = GR_SUCCESS;

  flint_rand_t state;
  flint_rand_init(state);

  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

  gr_mat_t A, G_a, H_a;
  gr_mat_t C, G_c, H_c;
  gr_mat_t I, G_i, H_i;
  gr_mat_t Check;
  gr_mat_init(A, n, m, ctx);
  gr_mat_init(I, n, m, ctx);
  gr_mat_init(C, n, m, ctx);
  gr_mat_init(Check, n, m, ctx);

  while (i < 20) {
    gr_mat_random_toeplitz(A, n, m, state, ctx);
    status |= gr_mat_one(I, ctx);

    status |= gr_mat_G_H(G_a, H_a, A, DISP_PLUS, ctx);
    status |= gr_mat_G_H(G_c, H_c, C, DISP_PLUS, ctx);
    status |= gr_mat_G_H(G_i, H_i, I, DISP_PLUS, ctx);

    status |= gr_mat_mul_generator(G_c, H_c, G_a, H_a, G_i, H_i, ctx);
    
    status |= gr_mat_generator_compress(G_c, H_c, G_c, H_c, ctx);
    
    status |= gr_mat_reconstruct_A(Check, G_c, H_c, DISP_PLUS, ctx);
    
    if (gr_mat_equal(Check, A, ctx) == T_FALSE){
        fprintf(stderr, "Test failed on the following matrices A:\n");
        gr_mat_print(A, ctx);
        gr_mat_print(Check, ctx);
        status = GR_TEST_FAIL;
        break;
    }
    i++;
  }

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
  if (strcmp("compression_then_multiply", argv[1]) == 0) {
    ok = test_compression_then_multiply();
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
