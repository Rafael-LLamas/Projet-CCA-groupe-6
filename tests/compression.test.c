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
  gr_mat_clear(Check, ctx);
  if (new_rank > old_rank) status = GR_TEST_FAIL;

  gr_mat_t D;
  gr_mat_init(D, n, m, ctx);
  status |= gr_mat_displacement(D, ref, DISP_PLUS, ctx);
  slong disp_rank;
  status |= gr_mat_rank(&disp_rank, D, ctx);
  gr_mat_clear(D, ctx);
  if (new_rank != disp_rank) status = GR_TEST_FAIL;
  return status;
}

int test_compress_toeplitz(int nb_iter) {
  flint_rand_t state;
  flint_rand_init(state);
  int status = GR_SUCCESS;
  int i = 0;
  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, n_randprime(state, 61, 1));

  while (i < nb_iter) {
    slong n = n_randint(state, 500);
    slong m = n_randint(state, 500);
    gr_mat_t A, G, H;
    gr_mat_init(A, n, m, ctx);
    status |= gr_mat_random_toeplitz(A, state, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("Error to create A");
      gr_mat_clear(A, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }
    status |= gr_mat_G_H(G, H, A, DISP_PLUS, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("Error to create G and H");
      gr_mat_clear(A, ctx);
      gr_mat_clear(G, ctx);
      gr_mat_clear(H, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }

    slong old_rank = gr_mat_ncols(G, ctx);
    status |= gr_mat_generator_compress(G, H, ctx);
    if (status != GR_SUCCESS) {
      flint_printf(" Error to Compress");
      gr_mat_clear(A, ctx);
      gr_mat_clear(G, ctx);
      gr_mat_clear(H, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }
    status |= check_compress_correctness(A, G, H, old_rank, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("failed verification compress");
      gr_mat_clear(A, ctx);
      gr_mat_clear(G, ctx);
      gr_mat_clear(H, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }

    gr_mat_clear(A, ctx);
    gr_mat_clear(G, ctx);
    gr_mat_clear(H, ctx);
    i++;
  }

  flint_rand_clear(state);
  gr_ctx_clear(ctx);
  return status;
}

int test_compress_after_addition(int nb_iter) {
  flint_rand_t state;
  flint_rand_init(state);
  int status = GR_SUCCESS;

  int i = 0;
  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, n_randprime(state, 60, 1));

  while (i < nb_iter) {
    slong n = n_randint(state, 500);
    slong m = n_randint(state, 500);
    gr_mat_t A, I, Ref, G_a, H_a, G_i, H_i, G_c, H_c;
    gr_mat_init(A, n, m, ctx);
    gr_mat_init(I, n, m, ctx);
    gr_mat_init(Ref, n, m, ctx);

    status |= gr_mat_random_toeplitz(A, state, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("Error to create A");
      gr_mat_clear(A, ctx);
      gr_mat_clear(I, ctx);
      gr_mat_clear(Ref, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }
    status |= gr_mat_one(I, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("Error to create I");
      gr_mat_clear(A, ctx);
      gr_mat_clear(I, ctx);
      gr_mat_clear(Ref, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }
    status |= gr_mat_add(Ref, A, I, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("Error to do the addition");
      gr_mat_clear(A, ctx);
      gr_mat_clear(I, ctx);
      gr_mat_clear(Ref, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }

    status |= gr_mat_G_H(G_a, H_a, A, DISP_PLUS, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("Error to create G_a and H_a");
      gr_mat_clear(A, ctx);
      gr_mat_clear(I, ctx);
      gr_mat_clear(Ref, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }
    status |= gr_mat_G_H(G_i, H_i, I, DISP_PLUS, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("Error to create G_i & H_i");
      gr_mat_clear(A, ctx);
      gr_mat_clear(I, ctx);
      gr_mat_clear(Ref, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(G_i, ctx);
      gr_mat_clear(H_i, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }
    status |= gr_mat_addition_generateur(G_a, H_a, G_i, H_i, G_c, H_c, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("Error addition the generators");
      gr_mat_clear(A, ctx);
      gr_mat_clear(I, ctx);
      gr_mat_clear(Ref, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(G_i, ctx);
      gr_mat_clear(H_i, ctx);
      gr_mat_clear(G_c, ctx);
      gr_mat_clear(H_c, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }

    slong old_rank = gr_mat_ncols(G_c, ctx);
    status |= gr_mat_generator_compress(G_c, H_c, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("Error to compress");
      gr_mat_clear(A, ctx);
      gr_mat_clear(I, ctx);
      gr_mat_clear(Ref, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(G_i, ctx);
      gr_mat_clear(H_i, ctx);
      gr_mat_clear(G_c, ctx);
      gr_mat_clear(H_c, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }
    status |= check_compress_correctness(Ref, G_c, H_c, old_rank, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("failed verify the compression");
      gr_mat_clear(A, ctx);
      gr_mat_clear(I, ctx);
      gr_mat_clear(Ref, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(G_i, ctx);
      gr_mat_clear(H_i, ctx);
      gr_mat_clear(G_c, ctx);
      gr_mat_clear(H_c, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }

    gr_mat_clear(A, ctx);
    gr_mat_clear(I, ctx);
    gr_mat_clear(Ref, ctx);
    gr_mat_clear(G_a, ctx);
    gr_mat_clear(H_a, ctx);
    gr_mat_clear(G_i, ctx);
    gr_mat_clear(H_i, ctx);
    gr_mat_clear(G_c, ctx);
    gr_mat_clear(H_c, ctx);
    i++;
  }

  flint_rand_clear(state);
  gr_ctx_clear(ctx);
  return status;
}

int test_compress_after_multiplication(int nb_iter) {
  flint_rand_t state;
  flint_rand_init(state);
  int status = GR_SUCCESS;

  int i = 0;

  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, n_randprime(state, 61, 1));

  while (i < nb_iter) {
    slong n = n_randint(state, 500);
    slong m = n_randint(state, 500);
    gr_mat_t A, I, G_a, H_a, G_i, H_i, G_c, H_c;
    gr_mat_init(A, n, m, ctx);
    gr_mat_init(I, m, n, ctx);

    status |= gr_mat_random_toeplitz(A, state, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("Error to create A");
      gr_mat_clear(A, ctx);
      gr_mat_clear(I, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }
    status |= gr_mat_one(I, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("Error to create I");
      gr_mat_clear(A, ctx);
      gr_mat_clear(I, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }

    status |= gr_mat_G_H(G_a, H_a, A, DISP_PLUS, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("Error to create G_a & H_a");
      gr_mat_clear(A, ctx);
      gr_mat_clear(I, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }
    status |= gr_mat_G_H(G_i, H_i, I, DISP_PLUS, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("Error to create G_i & H_i");
      gr_mat_clear(A, ctx);
      gr_mat_clear(I, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(G_i, ctx);
      gr_mat_clear(H_i, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }
    status |= gr_mat_mul_generator(G_c, H_c, G_a, H_a, G_i, H_i, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("Error to multiply the generators");
      gr_mat_clear(A, ctx);
      gr_mat_clear(I, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(G_i, ctx);
      gr_mat_clear(H_i, ctx);
      gr_mat_clear(G_c, ctx);
      gr_mat_clear(H_c, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }

    slong old_rank = gr_mat_ncols(G_c, ctx);
    status |= gr_mat_generator_compress(G_c, H_c, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("Error to crompress");
      gr_mat_clear(A, ctx);
      gr_mat_clear(I, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(G_i, ctx);
      gr_mat_clear(H_i, ctx);
      gr_mat_clear(G_c, ctx);
      gr_mat_clear(H_c, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }
    status |= check_compress_correctness(A, G_c, H_c, old_rank, ctx);
    if (status != GR_SUCCESS) {
      flint_printf("Failed to verify the compression");
      gr_mat_clear(A, ctx);
      gr_mat_clear(I, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(G_i, ctx);
      gr_mat_clear(H_i, ctx);
      gr_mat_clear(G_c, ctx);
      gr_mat_clear(H_c, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return status;
    }

    gr_mat_clear(A, ctx);
    gr_mat_clear(I, ctx);
    gr_mat_clear(G_a, ctx);
    gr_mat_clear(H_a, ctx);
    gr_mat_clear(G_i, ctx);
    gr_mat_clear(H_i, ctx);
    gr_mat_clear(G_c, ctx);
    gr_mat_clear(H_c, ctx);
    i++;
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
    ok = test_compress_toeplitz(30);
  } else if (strcmp("compress_after_addition", argv[1]) == 0) {
    ok = test_compress_after_addition(30);
  } else if (strcmp("compress_after_multiplication", argv[1]) == 0) {
    ok = test_compress_after_multiplication(30);
  } else {
    fprintf(stderr, "status: test \"%s\" not found!\n", argv[1]);
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
