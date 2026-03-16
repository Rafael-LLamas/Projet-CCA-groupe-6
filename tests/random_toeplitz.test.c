#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/ulong_extras.h"
#include "random_toeplitz.h"

int test_random_toepltiz() {
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  int i = 0;
  while (i < 10) {
    gr_ctx_t ctx;
    int error;
    gr_mat_t ran;

    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
    slong *rank = flint_malloc(sizeof(slong));
    gr_mat_init(ran, 5, 5, ctx);
    error = random_toeplitz(ran, 5, 5, state, ctx);
    if (error != 0) {
      gr_mat_clear(ran, ctx);
      gr_ctx_clear(ctx);
      flint_free(rank);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_rank(rank, ran, ctx);
    if (error != 0) {
      gr_mat_clear(ran, ctx);
      gr_ctx_clear(ctx);
      flint_free(rank);
      flint_rand_clear(state);
      return error;
    }
    flint_printf("----------------------------------------------\nResultat toeplitz = \n");
    gr_mat_print(ran, ctx);
    flint_printf("\n");
    flint_printf("rank = %wd \n", *rank);
    gr_mat_clear(ran, ctx);
    gr_ctx_clear(ctx);
    flint_free(rank);
    i++;
  }

  flint_rand_clear(state);
  return GR_SUCCESS;
}
int test_random_quasi_toeplitz() {
  int i = 0;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  while (i < 10) {
    gr_ctx_t ctx;
    int error;
    gr_mat_t ran, D;

    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
    slong *rank = flint_malloc(sizeof(slong));
    gr_mat_init(ran, 5, 5, ctx);
    error = rand_quasi_toeplitz(ran, 5, 5, 1, ctx);
    if (error != 0) {
      gr_mat_clear(ran, ctx);
      gr_ctx_clear(ctx);
      flint_free(rank);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_rank(rank, ran, ctx);
    if (error != 0) {
      flint_printf("%d\n", error);
      gr_mat_clear(ran, ctx);
      gr_ctx_clear(ctx);
      flint_free(rank);
      flint_rand_clear(state);
      return error;
    }
    flint_printf("----------------------------------------------\nResultat quasi 1 = \n");
    gr_mat_print(ran, ctx);
    flint_printf("\n");
    flint_printf("rank = %wd \n", *rank);
    gr_mat_clear(ran, ctx);
    gr_mat_init(ran, 5, 5, ctx);
    error = rand_quasi_toeplitz(ran, 5, 5, 4, ctx);
    if (error != 0) {
      gr_mat_clear(ran, ctx);
      gr_ctx_clear(ctx);
      flint_free(rank);
      flint_rand_clear(state);
      return error;
    }

    flint_printf("----------------------------------------------\nResultat quasi 4 = \n");
    gr_mat_print(ran, ctx);
    flint_printf("\n");
    gr_mat_clear(ran, ctx);
    gr_mat_init(ran, 5, 5, ctx);
    error = rand_quasi_toeplitz(ran, 5, 5, 2, ctx);
    if (error != 0) {
      gr_mat_clear(ran, ctx);
      gr_ctx_clear(ctx);
      flint_free(rank);
      flint_rand_clear(state);
      return error;
    }
    flint_printf("----------------------------------------------\nResultat quasi 2 = \n");
    gr_mat_print(ran, ctx);
    flint_printf("\n");
    gr_mat_clear(ran, ctx);
    gr_mat_init(ran, 5, 5, ctx);
    error = rand_quasi_toeplitz(ran, 5, 5, 3, ctx);
    if (error != 0) {
      gr_mat_clear(ran, ctx);
      gr_ctx_clear(ctx);
      flint_free(rank);
      flint_rand_clear(state);
      return error;
    }

    flint_printf("----------------------------------------------\nResultat quasi 3 = \n");
    gr_mat_print(ran, ctx);
    flint_printf("\n");
    gr_mat_clear(ran, ctx);
    gr_mat_init(ran, 15, 15, ctx);
    error = rand_quasi_toeplitz(ran, 15, 15, 0, ctx);
    if (error != 0) {
      gr_mat_clear(ran, ctx);
      gr_ctx_clear(ctx);
      flint_free(rank);
      flint_rand_clear(state);
      return error;
    }
    gr_mat_init(D, 15, 15, ctx);
    error = gr_mat_displacement_square_safe(D, ran, ctx);
    if (error != 0) {
      gr_mat_clear(ran, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_free(rank);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_rank(rank, D, ctx);
    if (error != 0) {
      gr_mat_clear(ran, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_free(rank);
      flint_rand_clear(state);
      return error;
    }
    flint_printf("----------------------------------------------\nResultat quasi 0 = \n");
    gr_mat_print(ran, ctx);
    flint_printf("\n\n\n");
    flint_printf(" deplacement =  ");
    gr_mat_print(D, ctx);
    flint_printf("\n");
    flint_printf("rank de déplcement = %wd \n", *rank);

    gr_mat_clear(D, ctx);
    gr_mat_clear(ran, ctx);
    gr_ctx_clear(ctx);
    flint_free(rank);
    i++;
  }

  flint_rand_clear(state);
  return GR_SUCCESS;
}

void usage(char *argv[]) {
  fprintf(stderr, "Usage: %s <test_name>\n", argv[0]);
  fprintf(stderr, "Available tests:\n");
  fprintf(stderr, "  - random_toeplitz\n");
  fprintf(stderr, "  - random_quasi_toeplitz\n");
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    usage(argv);
    return GR_UNABLE;
  }
  // start test
  fprintf(stderr, "=> Start test \"%s\"\n", argv[1]);
  int ok = GR_SUCCESS;
  if (strcmp("random_toeplitz", argv[1]) == 0) {
    ok = test_random_toepltiz();
  } else if (strcmp("random_quasi_toeplitz", argv[1]) == 0) {
    ok = test_random_quasi_toeplitz();
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
