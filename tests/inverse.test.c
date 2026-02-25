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
#include "toeplitz_inverse.h"

int test_inverse_base_case() {
  int res = GR_SUCCESS;
  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, 17);
  slong n = 1;
  slong rank = 2;
  gr_mat_t G_A, H_A, G_D, H_D;
  gr_mat_init(G_A, n, rank, ctx);
  gr_mat_init(H_A, n, rank, ctx);
  gr_mat_init(G_D, n, rank, ctx);
  gr_mat_init(H_D, n, rank, ctx);
  gr_set_si(gr_mat_entry_ptr(G_A, 0, 0, ctx), 3, ctx);
  gr_set_si(gr_mat_entry_ptr(G_A, 0, 1, ctx), -1, ctx);
  gr_set_si(gr_mat_entry_ptr(H_A, 0, 0, ctx), 2, ctx);
  gr_set_si(gr_mat_entry_ptr(H_A, 0, 1, ctx), 1, ctx);
  flint_printf("Testing Strassen Inverse Base Case (n=1) over nmod 17\n");
  flint_printf("---------------------------------------------------\n");
  flint_printf("Input G_A: ");
  gr_mat_print(G_A, ctx);
  flint_printf("\nInput H_A: ");
  gr_mat_print(H_A, ctx);
  flint_printf("\n\n");
  int out = gr_toeplitz_inverse(G_D, H_D, G_A, H_A, ctx);
  if (out == GR_SUCCESS) {
    flint_printf("[SUCCESS] Matrix inverted.\n");
    flint_printf("Expected scalar: 3*2 + (-1)*1 = 5\n");
    flint_printf("Expected inverse of 5 mod 17: 7\n\n");

    flint_printf("Output G_D: ");
    gr_mat_print(G_D, ctx); // Should be [7, 0]
    flint_printf("\nOutput H_D: ");
    gr_mat_print(H_D, ctx); // Should be [1, 0]
    flint_printf("\n");
  } else {
    flint_printf("[ERROR] Failed to invert matrix.\n");
    res = GR_TEST_FAIL;
  }
  gr_mat_clear(G_A, ctx);
  gr_mat_clear(H_A, ctx);
  gr_mat_clear(G_D, ctx);
  gr_mat_clear(H_D, ctx);
  gr_ctx_clear(ctx);
  return res;
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
  if (strcmp("inverse_base_case", argv[1]) == 0) {
    ok = test_inverse_base_case();
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