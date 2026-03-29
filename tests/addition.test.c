#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "addition.h"
#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/ulong_extras.h"
#include "matrix_aux.h"
#include "random_toeplitz.h"

int test_addition_generateurs_toeplitz_square() {
  int i = 0;
  flint_rand_t state;
  flint_rand_init(state);
  while (i < 10) {
    gr_mat_t A1, B1, C1;
    gr_ctx_t ctx;
    flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
    int N = 5;
    flint_printf("Reference (FLINT):\n");
    gr_mat_init(A1, N, N, ctx);
    gr_mat_init(B1, N, N, ctx);
    gr_mat_init(C1, N, N, ctx);
    FLINT_CHECK(gr_mat_random_toeplitz(A1, state, ctx));
    FLINT_CHECK(gr_mat_random_toeplitz(B1, state, ctx));
    FLINT_CHECK(gr_mat_add(C1, A1, B1, ctx));
    flint_printf("A1 = \n");
    gr_mat_print(A1, ctx);
    flint_printf("\n");
    flint_printf("B1 = \n");
    gr_mat_print(B1, ctx);
    flint_printf("\n");
    flint_printf("C1 = \n");
    gr_mat_print(C1, ctx);
    flint_printf("\n");
    flint_printf("With Generators:\n");
    gr_mat_t A1_G, A1_H, B1_G, B1_H, C2_G, C2_H, C2;
    gr_mat_init(C2, N, N, ctx);

    FLINT_CHECK(gr_mat_G_H(A1_G, A1_H, A1, DISP_PLUS, ctx));
    FLINT_CHECK(gr_mat_G_H(B1_G, B1_H, B1, DISP_PLUS, ctx));
    FLINT_CHECK(gr_mat_addition_generateur(A1_G, A1_H, B1_G, B1_H, C2_G, C2_H, ctx));
    FLINT_CHECK(gr_mat_reconstruct_A(C2, C2_G, C2_H, DISP_PLUS, ctx));

    flint_printf("A1_G = \n");
    gr_mat_print(A1_G, ctx);
    flint_printf("\n");
    flint_printf("A1_H = \n");
    gr_mat_print(A1_H, ctx);
    flint_printf("\n");
    flint_printf("B1 = \n");
    gr_mat_print(B1_G, ctx);
    flint_printf("\n");
    flint_printf("B1_H = \n");
    gr_mat_print(B1_H, ctx);
    flint_printf("\n");
    flint_printf("C2 = \n");
    gr_mat_print(C2, ctx);
    flint_printf("\n");

    FLINT_CHECK(gr_mat_equal(C1, C2, ctx));

    gr_mat_clear(A1, ctx);
    gr_mat_clear(B1, ctx);
    gr_mat_clear(C1, ctx);
    gr_mat_clear(A1_G, ctx);
    gr_mat_clear(A1_H, ctx);
    gr_mat_clear(B1_G, ctx);
    gr_mat_clear(B1_H, ctx);
    gr_mat_clear(C2_G, ctx);
    gr_mat_clear(C2_H, ctx);
    gr_mat_clear(C2, ctx);
    gr_ctx_clear(ctx);
    i++;
  }
  flint_rand_clear(state);
  return GR_SUCCESS;
}

void usage(char *argv[]) {
  fprintf(stderr, "Usage: %s <test_name>\n", argv[0]);
  fprintf(stderr, "Available tests:\n");
  fprintf(stderr, "  - addition_generators_toeplitz\n");
}
int main(int argc, char *argv[]) {
  if (argc < 2) {
    usage(argv);
    return GR_UNABLE;
  }

  fprintf(stderr, "=> Start test \"%s\"\n", argv[1]);
  int ok = GR_SUCCESS;

  if (strcmp("addition_generators_toeplitz", argv[1]) == 0) {
    ok = test_addition_generateurs_toeplitz_square();
  } else {
    fprintf(stderr, "Error: test \"%s\" not found!\n", argv[1]);
    exit(EXIT_FAILURE);
  }

  if (ok == GR_SUCCESS) {
    fprintf(stderr, "Test \"%s\" finished: SUCCESS\n", argv[1]);
    return EXIT_SUCCESS;
  } else {
    fprintf(stderr, "Test \"%s\" finished: FAILURE\n", argv[1]);
    return EXIT_FAILURE;
  }
}
