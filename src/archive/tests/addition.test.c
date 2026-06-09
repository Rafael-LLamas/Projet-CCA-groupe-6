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

int test_addition_generateurs_toeplitz() {
  int i = 0;
  int error = GR_SUCCESS;
  flint_rand_t state;
  flint_rand_init(state);
  while (i < 50) {
    gr_mat_t A, B, C;
    gr_ctx_t ctx;
    flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
    slong n = n_randint(state, 100);
    slong m = n_randint(state, 100);

    gr_mat_init(A, n, m, ctx);
    gr_mat_init(B, n, m, ctx);
    gr_mat_init(C, n, m, ctx);
    error = gr_mat_random_toeplitz(A, state, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("A is not toeplitz, fix it");
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_random_toeplitz(B, state, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("B is not toeplitz, fix it");
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_add(C, A, B, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Bad Matrix ?");
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    gr_mat_t A_G, A_H, B_G, B_H, D_G, D_H, D;
    gr_mat_init(D, n, m, ctx);

    error = gr_mat_G_H(A_G, A_H, A, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf(" Cant get G and H of A !!");
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_mat_clear(A_G, ctx);
      gr_mat_clear(A_H, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_G_H(B_G, B_H, B, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf(" Cant get G and H of B");
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_mat_clear(A_G, ctx);
      gr_mat_clear(A_H, ctx);
      gr_mat_clear(B_G, ctx);
      gr_mat_clear(B_H, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_addition_generateur(A_G, A_H, B_G, B_H, D_G, D_H, ctx);
    if (error != GR_SUCCESS) {
      flint_printf(" How !? that just a concatenation !!");
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_mat_clear(A_G, ctx);
      gr_mat_clear(A_H, ctx);
      gr_mat_clear(B_G, ctx);
      gr_mat_clear(B_H, ctx);
      gr_mat_clear(D_G, ctx);
      gr_mat_clear(D_H, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_reconstruct_A(D, D_G, D_H, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf(" Cant reconstruct, fix it");
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_mat_clear(A_G, ctx);
      gr_mat_clear(A_H, ctx);
      gr_mat_clear(B_G, ctx);
      gr_mat_clear(B_H, ctx);
      gr_mat_clear(D_G, ctx);
      gr_mat_clear(D_H, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    if (gr_mat_equal(C, D, ctx) != T_TRUE) {
      flint_printf(" C NOT EQUAL D !!!\n");
      flint_printf("A = \n");
      gr_mat_print(A, ctx);
      flint_printf("\n");
      flint_printf("B = \n");
      gr_mat_print(B, ctx);
      flint_printf("\n");
      flint_printf("C = \n");
      gr_mat_print(C, ctx);
      flint_printf("\n");
      flint_printf("D = \n");
      gr_mat_print(D, ctx);
      flint_printf("\n");
      flint_printf("A_G = \n");
      gr_mat_print(A_G, ctx);
      flint_printf("\n");
      flint_printf("A_H = \n");
      gr_mat_print(A_H, ctx);
      flint_printf("\n");
      flint_printf("B_G = \n");
      gr_mat_print(B_G, ctx);
      flint_printf("\n");
      flint_printf("B_H = \n");
      gr_mat_print(B_H, ctx);
      flint_printf("\n");
      flint_printf("D_G = \n");
      gr_mat_print(D_G, ctx);
      flint_printf("\n");
      flint_printf("D_H = \n");
      gr_mat_print(D_H, ctx);
      flint_printf("\n");
      i = 9998;
    }

    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);
    gr_mat_clear(A_G, ctx);
    gr_mat_clear(A_H, ctx);
    gr_mat_clear(B_G, ctx);
    gr_mat_clear(B_H, ctx);
    gr_mat_clear(D_G, ctx);
    gr_mat_clear(D_H, ctx);
    gr_mat_clear(D, ctx);
    gr_ctx_clear(ctx);
    i++;
  }
  flint_rand_clear(state);
  return error;
}
int test_addition_generateurs_quasi_toeplitz() {
  int i = 0;
  int error = GR_SUCCESS;
  flint_rand_t state;
  flint_rand_init(state);
  while (i < 50) {
    gr_mat_t A, B, C;
    gr_ctx_t ctx;
    flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
    slong n = n_randint(state, 100);
    slong m = n_randint(state, 100);
    slong r = n_randint(state, 5);

    gr_mat_init(A, n, m, ctx);
    gr_mat_init(B, n, m, ctx);
    gr_mat_init(C, n, m, ctx);

    error = gr_mat_quasi_toeplitz_rank(A, r, state, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("A is not quasi toeplitz, fix it");
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_quasi_toeplitz_rank(B, r, state, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("B is not quasi toeplitz, fix it");
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_add(C, A, B, ctx);
    if (error != GR_SUCCESS) {
      flint_printf("Bad Matrix ?");
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    gr_mat_t A_G, A_H, B_G, B_H, D_G, D_H, D;
    gr_mat_init(D, n, m, ctx);

    error = gr_mat_G_H(A_G, A_H, A, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf(" Cant get G and H of A !!");
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_mat_clear(A_G, ctx);
      gr_mat_clear(A_H, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_G_H(B_G, B_H, B, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf(" Cant get G and H of B");
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_mat_clear(A_G, ctx);
      gr_mat_clear(A_H, ctx);
      gr_mat_clear(B_G, ctx);
      gr_mat_clear(B_H, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_addition_generateur(A_G, A_H, B_G, B_H, D_G, D_H, ctx);
    if (error != GR_SUCCESS) {
      flint_printf(" How !? that just a concatenation !!\n");
      flint_printf("rows of G_A = %d\n", gr_mat_nrows(A_G, ctx));
      flint_printf("rows of H_A = %d\n", gr_mat_nrows(A_H, ctx));
      flint_printf("rows of G_B = %d\n", gr_mat_nrows(B_G, ctx));
      flint_printf("rows of H_B = %d\n", gr_mat_nrows(B_H, ctx));
      flint_printf("cols of G_A = %d\n", gr_mat_ncols(A_G, ctx));
      flint_printf("cols of H_A = %d\n", gr_mat_ncols(A_H, ctx));
      flint_printf("cols of G_B = %d\n", gr_mat_ncols(B_G, ctx));
      flint_printf("cols of H_B = %d\n", gr_mat_ncols(B_H, ctx));
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_mat_clear(A_G, ctx);
      gr_mat_clear(A_H, ctx);
      gr_mat_clear(B_G, ctx);
      gr_mat_clear(B_H, ctx);
      gr_mat_clear(D_G, ctx);
      gr_mat_clear(D_H, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    error = gr_mat_reconstruct_A(D, D_G, D_H, DISP_PLUS, ctx);
    if (error != GR_SUCCESS) {
      flint_printf(" Cant reconstruct, fix it\n");
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_mat_clear(A_G, ctx);
      gr_mat_clear(A_H, ctx);
      gr_mat_clear(B_G, ctx);
      gr_mat_clear(B_H, ctx);
      gr_mat_clear(D_G, ctx);
      gr_mat_clear(D_H, ctx);
      gr_mat_clear(D, ctx);
      gr_ctx_clear(ctx);
      flint_rand_clear(state);
      return error;
    }
    if (gr_mat_equal(C, D, ctx) != T_TRUE) {
      flint_printf(" C NOT EQUAL D !!!\n");
      flint_printf("A = \n");
      gr_mat_print(A, ctx);
      flint_printf("\n");
      flint_printf("B = \n");
      gr_mat_print(B, ctx);
      flint_printf("\n");
      flint_printf("C = \n");
      gr_mat_print(C, ctx);
      flint_printf("\n");
      flint_printf("D = \n");
      gr_mat_print(D, ctx);
      flint_printf("\n");
      flint_printf("A_G = \n");
      gr_mat_print(A_G, ctx);
      flint_printf("\n");
      flint_printf("A_H = \n");
      gr_mat_print(A_H, ctx);
      flint_printf("\n");
      flint_printf("B_G = \n");
      gr_mat_print(B_G, ctx);
      flint_printf("\n");
      flint_printf("B_H = \n");
      gr_mat_print(B_H, ctx);
      flint_printf("\n");
      flint_printf("D_G = \n");
      gr_mat_print(D_G, ctx);
      flint_printf("\n");
      flint_printf("D_H = \n");
      gr_mat_print(D_H, ctx);
      flint_printf("\n");
      i = 9998;
    }

    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);
    gr_mat_clear(A_G, ctx);
    gr_mat_clear(A_H, ctx);
    gr_mat_clear(B_G, ctx);
    gr_mat_clear(B_H, ctx);
    gr_mat_clear(D_G, ctx);
    gr_mat_clear(D_H, ctx);
    gr_mat_clear(D, ctx);
    gr_ctx_clear(ctx);
    i++;
  }
  flint_rand_clear(state);
  return error;
}

void usage(char *argv[]) {
  fprintf(stderr, "Usage: %s <test_name>\n", argv[0]);
  fprintf(stderr, "Available tests:\n");
  fprintf(stderr, "  - addition_generators_toeplitz\n");
  fprintf(stderr, "  - addition_generators_quasi_toeplitz\n");
}
int main(int argc, char *argv[]) {
  if (argc < 2) {
    usage(argv);
    return GR_UNABLE;
  }

  fprintf(stderr, "=> Start test \"%s\"\n", argv[1]);
  int ok = GR_SUCCESS;

  if (strcmp("addition_generators_toeplitz", argv[1]) == 0) {
    ok = test_addition_generateurs_toeplitz();
  } else if (strcmp("addition_generators_quasi_toeplitz", argv[1]) == 0) {
    ok = test_addition_generateurs_quasi_toeplitz();
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
