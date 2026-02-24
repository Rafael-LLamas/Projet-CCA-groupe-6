#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "flint/ulong_extras.h"
#include "matrix_aux.h"
#include "multiplication.h"
#include "random_toeplitz.h"

int test_multiplication_generateurs() {
  int error;
  gr_mat_t A, B, C, G_a, H_a, G_b, H_b, G_c, H_c;
  gr_ctx_t ctx;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
  gr_mat_init(A, 5, 5, ctx);
  gr_mat_init(B, 5, 5, ctx);
  error = random_toeplitz(A, 5, 5, state, ctx);
  if (error != 0) {
    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_ctx_clear(ctx);
    flint_rand_clear(state);
    return error;
  }
  error = random_toeplitz(B, 5, 5, state, ctx);
  if (error != 0) {

    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_ctx_clear(ctx);
    flint_rand_clear(state);
    return error;
  }
  flint_printf("A = \n");
  gr_mat_print(A, ctx);
  flint_printf("\n\n");
  flint_printf("B = \n");
  gr_mat_print(B, ctx);
  flint_printf("\n\n");
  gr_mat_init(C, gr_mat_nrows(A, ctx), gr_mat_ncols(B, ctx), ctx);
  error = gr_mat_G_H(G_a, H_a, A, ctx);
  if (error != 0) {
    gr_mat_clear(C, ctx);
    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(G_a, ctx);
    gr_mat_clear(H_a, ctx);
    gr_ctx_clear(ctx);
    flint_rand_clear(state);
    return error;
  }
  error = gr_mat_G_H(G_b, H_b, B, ctx);
  if (error != 0) {
    gr_mat_clear(C, ctx);
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
  if (error != 0) {
    gr_mat_clear(C, ctx);
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
  flint_printf(":------------------------Matrices C faite sans les générateurs---------------------------------:\n");
  flint_printf("Matrice C = \n");
  gr_mat_print(C, ctx);
  flint_printf("\n");
  flint_printf(":------------------------Matrices C faite avec les générateurs---------------------------------:\n");
  error = gr_multiplication_generateur_deplacement_fast(G_c, H_c, G_a, H_a, G_b, H_b, ctx);
  if (error != 0) {
    gr_mat_clear(C, ctx);
    gr_mat_clear(A, ctx);
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
  error = gr_mat_reconstruct_A_safe(C, G_c, H_c, ctx);
  if (error != 0) {
    gr_mat_clear(C, ctx);
    gr_mat_clear(A, ctx);
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
  flint_printf("Matrice C = \n");
  gr_mat_print(C, ctx);
  flint_printf("\n");

  gr_mat_clear(C, ctx);
  gr_mat_clear(A, ctx);
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
void usage(char *argv[]) {
  fprintf(stderr, "Usage: %s <test_name>\n", argv[0]);
  fprintf(stderr, "Available tests:\n");
  fprintf(stderr, "  - multiplication_generators\n");
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    usage(argv);
    return GR_UNABLE;
  }
  // start test
  fprintf(stderr, "=> Start test \"%s\"\n", argv[1]);
  int ok = GR_SUCCESS;
  if (strcmp("multiplication_generators", argv[1]) == 0) {
    ok = test_multiplication_generateurs();
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