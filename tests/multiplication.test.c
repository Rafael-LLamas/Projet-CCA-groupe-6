#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "matrix_aux.h"
#include "multiplication.h"
#include "random_toeplitz.h"

int test_multiplication_generateurs() {
  int error;
  gr_mat_t A, B, C, G_a, H_a, G_b, H_b, G_c, H_c;
  gr_ctx_t ctx;
  flint_rand_t state;
  gr_ctx_init_nmod(ctx, 47);
  flint_rand_init(state);
  error = random_toeplitz(A, gr_mat_nrows(A, ctx), gr_mat_ncols(A, ctx), state, ctx);
  if (error != 0) { return error; }
  error = random_toeplitz(B, gr_mat_nrows(B, ctx), gr_mat_ncols(B, ctx), state, ctx);
  if (error != 0) { return error; }
  gr_mat_init(C, gr_mat_nrows(A, ctx), gr_mat_ncols(B, ctx), ctx);
  error = gr_mat_G_H(G_a, H_a, A, ctx);
  if (error != 0) { return error; }
  error = gr_mat_G_H(G_b, H_b, B, ctx);
  if (error != 0) { return error; }
  error = gr_multiplication_toeplitz(C, A, B, ctx);
  if (error != 0) { return error; }
  flint_printf(":------------------------Matrices C faite sans les générateurs---------------------------------:");
  flint_printf("Matrice C = \n");
  gr_mat_print(C, ctx);
  flint_printf("\n");
  gr_mat_G_H(C, G_c, H_c, ctx);
  flint_printf("Matrice G_c = \n");
  gr_mat_print(G_c, ctx);
  flint_printf("\n");
  flint_printf("Matrice H_c = \n");
  gr_mat_print(H_c, ctx);
  flint_printf("\n");
  gr_multiplication_generateur_deplacement(G_c, H_c, G_a, H_a, G_b, H_b, ctx);
  flint_printf(":------------------------Matrices C faite avec les générateurs---------------------------------:");
  flint_printf("Matrice G_c = \n");
  gr_mat_print(G_c, ctx);
  flint_printf("\n");
  flint_printf("Matrice H_c = \n");
  gr_mat_print(H_c, ctx);
  flint_printf("\n");
  gr_mat_reconstruct_A_safe(C, G_c, H_c, ctx);
  flint_printf("Matrice C = \n");
  gr_mat_print(C, ctx);
  flint_printf("\n");
  return error;
}

int main(int argc, char *argv[]) {
  if (argc == 1) { return GR_UNABLE; }
  // start test
  fprintf(stderr, "=> Start test \"%s\"\n", argv[1]);
  int ok = GR_SUCCESS;
  if (strcmp("multiplication_generators", argv[1]) == 0) {
    ok = test_multiplication_generateurs();
  } else if (strcmp("this is your DIY project rafael hamas", argv[1]) == 0) {
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