#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "addition.h"
#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "matrix_aux.h"
#include "random_toeplitz.h"

int test_addition_generateurs_toeplitz_square() {
  gr_mat_t A1, B1, C1;
  gr_ctx_t ctx;
  flint_rand_t state;
  gr_ctx_init_nmod(ctx, GNMOD);
  flint_rand_init(state);
  int N = 5;
  flint_printf("Reference (FLINT):\n");
  gr_mat_init(A1, N, N, ctx);
  gr_mat_init(B1, N, N, ctx);
  gr_mat_init(C1, N, N, ctx);
  flint_rand_init(state);
  FLINT_CHECK(random_toeplitz(A1, gr_mat_nrows(A1, ctx), gr_mat_ncols(A1, ctx), state, ctx));
  FLINT_CHECK(random_toeplitz(B1, gr_mat_nrows(B1, ctx), gr_mat_ncols(B1, ctx), state, ctx));
  FLINT_CHECK(gr_mat_add(C1, A1, B1, ctx));
  gr_mat_print(C1, ctx);
  flint_printf("\nWith Generators:\n");
  gr_mat_t A1_G, A1_H, B1_G, B1_H, C2_G, C2_H, C2;
  gr_mat_init(C2, N, N, ctx);

  FLINT_CHECK(gr_mat_G_H(A1_G, A1_H, A1, ctx));
  FLINT_CHECK(gr_mat_G_H(B1_G, B1_H, B1, ctx));
  gr_mat_init(C2_G, N, 4, ctx); // moddif ici pour avoir un code dynamique mais normalement le G d'une toeplitz est de
                                // taille Nx2 donc le G_res = Nx4
  gr_mat_init(C2_H, N, 4, ctx);
  FLINT_CHECK(gr_mat_addition_generateur(A1_G, A1_H, B1_G, B1_H, C2_G, C2_H, ctx));
  FLINT_CHECK(gr_mat_reconstruct_A_safe(C2, C2_G, C2_H, ctx));
  gr_mat_print(C2, ctx);
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
  flint_rand_clear(state);
  return GR_SUCCESS;
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <test_name>\n", argv[0]);
    return GR_UNABLE;
  }

  fprintf(stderr, "=> Start test \"%s\"\n", argv[1]);
  int ok = GR_SUCCESS;

  if (strcmp("addition_generators_toeplitz", argv[1]) == 0) {
    ok = test_addition_generateurs_toeplitz_square();
  } else if (strcmp("placeholder", argv[1]) == 0) {
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