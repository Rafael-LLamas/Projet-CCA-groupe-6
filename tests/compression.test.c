#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "addition.h"
#include "compression.h"
#include "displacement_matrices.h"
#include "multiplication.h"
#include "random_toeplitz.h"

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"
#include "flint/ulong_extras.h"

// aux
int check_compress_correctness(gr_mat_t ref, gr_mat_t G_c, gr_mat_t H_c, slong old_rank, gr_ctx_t ctx) {
  int status = GR_SUCCESS;
  slong n = gr_mat_nrows(ref, ctx);
  slong m = gr_mat_ncols(ref, ctx);
  slong new_rank = gr_mat_ncols(G_c, ctx);
  gr_mat_t Check;
  gr_mat_init(Check, n, m, ctx);
  status |= gr_mat_reconstruct_A(Check, G_c, H_c, DISP_PLUS, ctx);
  if (gr_mat_equal(Check, ref, ctx) == T_FALSE) status = GR_TEST_FAIL;
  if (new_rank > old_rank) status = GR_TEST_FAIL;
  gr_mat_clear(Check, ctx);
  return status;
}

int test_compress_toeplitz() {
  int status = GR_SUCCESS;
  int n = 20;

  flint_rand_t state;
  flint_rand_init(state);

  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, n_randprime(state, 61, 1));

  for (int i = 0; i < 30; i++) {
    gr_mat_t A, G, H;
    gr_mat_init(A, n, n, ctx);
    status |= gr_mat_random_toeplitz(A, n, n, state, ctx);
    status |= gr_mat_G_H(G, H, A, DISP_PLUS, ctx);

    slong old_rank = gr_mat_ncols(G, ctx);
    status |= gr_mat_generator_compress(G, H, ctx);
    status |= check_compress_correctness(A, G, H, old_rank, ctx);

    if (gr_mat_ncols(G, ctx) > 2) status = GR_TEST_FAIL;

    gr_mat_clear(A, ctx);
    gr_mat_clear(G, ctx);
    gr_mat_clear(H, ctx);
  }

  flint_rand_clear(state);
  gr_ctx_clear(ctx);
  return status;
}

int test_compress_after_addition() {
  int status = GR_SUCCESS;
  int n = 20;

  flint_rand_t state;
  flint_rand_init(state);

  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, n_randprime(state, 61, 1));

  for (int i = 0; i < 30; i++) {
    gr_mat_t A, I, Ref, G_a, H_a, G_i, H_i, G_c, H_c;
    gr_mat_init(A, n, n, ctx);
    gr_mat_init(I, n, n, ctx);
    gr_mat_init(Ref, n, n, ctx);

    status |= gr_mat_random_toeplitz(A, n, n, state, ctx);
    status |= gr_mat_one(I, ctx);
    status |= gr_mat_add(Ref, A, I, ctx);

    status |= gr_mat_G_H(G_a, H_a, A, DISP_PLUS, ctx);
    status |= gr_mat_G_H(G_i, H_i, I, DISP_PLUS, ctx);
    status |= gr_mat_addition_generateur(G_a, H_a, G_i, H_i, G_c, H_c, ctx);

    slong old_rank = gr_mat_ncols(G_c, ctx);
    status |= gr_mat_generator_compress(G_c, H_c, ctx);
    status |= check_compress_correctness(Ref, G_c, H_c, old_rank, ctx);

    gr_mat_clear(A, ctx);
    gr_mat_clear(I, ctx);
    gr_mat_clear(Ref, ctx);
    gr_mat_clear(G_a, ctx);
    gr_mat_clear(H_a, ctx);
    gr_mat_clear(G_i, ctx);
    gr_mat_clear(H_i, ctx);
    gr_mat_clear(G_c, ctx);
    gr_mat_clear(H_c, ctx);
  }

  flint_rand_clear(state);
  gr_ctx_clear(ctx);
  return status;
}

int test_compress_after_multiplication() {
  int status = GR_SUCCESS;
  int n = 20;

  flint_rand_t state;
  flint_rand_init(state);
  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, n_randprime(state, 61, 1));

  for (int i = 0; i < 30; i++) {
    gr_mat_t A, I, G_a, H_a, G_i, H_i, G_c, H_c;
    gr_mat_init(A, n, n, ctx);
    gr_mat_init(I, n, n, ctx);

    status |= gr_mat_random_toeplitz(A, n, n, state, ctx);
    status |= gr_mat_one(I, ctx);

    status |= gr_mat_G_H(G_a, H_a, A, DISP_PLUS, ctx);
    status |= gr_mat_G_H(G_i, H_i, I, DISP_PLUS, ctx);
    status |= gr_mat_mul_generator(G_c, H_c, G_a, H_a, G_i, H_i, ctx);

    slong old_rank = gr_mat_ncols(G_c, ctx);
    status |= gr_mat_generator_compress(G_c, H_c, ctx);
    status |= check_compress_correctness(A, G_c, H_c, old_rank, ctx);

    gr_mat_clear(A, ctx);
    gr_mat_clear(I, ctx);
    gr_mat_clear(G_a, ctx);
    gr_mat_clear(H_a, ctx);
    gr_mat_clear(G_i, ctx);
    gr_mat_clear(H_i, ctx);
    gr_mat_clear(G_c, ctx);
    gr_mat_clear(H_c, ctx);
  }

  flint_rand_clear(state);
  gr_ctx_clear(ctx);
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
  if (strcmp("compress_toeplitz", argv[1]) == 0) {
    ok = test_compress_toeplitz();
  } else if (strcmp("compress_after_addition", argv[1]) == 0) {
    ok = test_compress_after_addition();
  } else if (strcmp("compress_after_multiplication", argv[1]) == 0) {
    ok = test_compress_after_multiplication();
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
