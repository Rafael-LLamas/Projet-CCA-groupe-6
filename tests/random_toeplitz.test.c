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
    int error = GR_SUCCESS;
    gr_mat_t ran, D;
    slong n = n_randint(state, 1000);
    slong m = n_randint(state, 1000);
    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
    gr_mat_init(ran, n, m, ctx);
    gr_mat_init(D, n, m, ctx);
    error = gr_mat_random_toeplitz(ran, state, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to create a toeplitz, FIX IT FAST !\n");
      gr_mat_clear(ran, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_displacement(D, ran, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to do the displacement, FIX IT FAST !\n");
      gr_mat_clear(ran, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    slong r;
    error = gr_mat_rank(&r, D, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to get the rank, dunno why FLINT WHY FLINT\n");
      gr_mat_clear(ran, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    if (r != 2) {
      flint_printf("BAD RANK !!!!!\n");
      flint_printf("Matrice RAN = \n");
      gr_mat_print(ran, ctx);
      flint_printf("\n");
      flint_printf("Matrice D = \n");
      gr_mat_print(D, ctx);
      flint_printf("\n");
      flint_printf("rank displacement = %d\n", r);
      i = 9998;
      error = 1;
    }
    gr_mat_clear(D, ctx);
    gr_mat_clear(ran, ctx);
    gr_ctx_clear(ctx);
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
    int error = GR_SUCCESS;
    gr_mat_t ran, D;
    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
    slong n = n_randint(state, 1000);
    slong m = n_randint(state, 1000);
    slong rank = n_randint(state, 5);
    gr_mat_init(ran, n, m, ctx);
    gr_mat_init(D, n, m, ctx);
    error = gr_mat_quasi_toeplitz_rank(ran, rank, state, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to create a quasi toeplitz, FIX IT FAST !\n");
      gr_mat_clear(ran, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_displacement(D, ran, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to do the displacement, FIX IT FAST !\n");
      gr_mat_clear(ran, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    slong r;
    error = gr_mat_rank(&r, D, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Failed to get the rank, dunno why FLINT WHY FLINT\n");
      gr_mat_clear(ran, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    if (rank > 1 && r != rank) {
      flint_printf("BAD RANK !!!!!\n");
      flint_printf("Matrice RAN = \n");
      gr_mat_print(ran, ctx);
      flint_printf("\n");
      flint_printf("Matrice D = \n");
      gr_mat_print(D, ctx);
      flint_printf("\n");
      flint_printf("starting rank = %d\n", rank);
      flint_printf("rank displacement = %d\n", r);
      i = 9998;
      error = 1;
    }
    gr_mat_clear(D, ctx);
    gr_mat_clear(ran, ctx);
    gr_ctx_clear(ctx);
    i++;
  }

  flint_rand_clear(state);
  return GR_SUCCESS;
}

void usage(char *argv[]) {
  fprintf(stderr, "Usage: %s <test_name>\n", argv[0]);
  fprintf(stderr, "Available tests:\n");
  fprintf(stderr, "  - gr_mat_random_toeplitz\n");
  fprintf(stderr, "  - gr_mat_random_quasi_toepitz\n");
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    usage(argv);
    return GR_UNABLE;
  }
  // start test
  fprintf(stderr, "=> Start test \"%s\"\n", argv[1]);
  int ok = GR_SUCCESS;
  if (strcmp("gr_mat_random_toeplitz", argv[1]) == 0) {
    ok = test_random_toepltiz();
  } else if (strcmp("gr_mat_random_quasi_toeplitz", argv[1]) == 0) {
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
