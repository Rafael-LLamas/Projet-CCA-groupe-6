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

int test_inverse_2x2() {
  int res = GR_SUCCESS;
  gr_ctx_t ctx;
  ulong p = 97;
  gr_ctx_init_nmod(ctx, p);

  flint_printf("Testing Strassen Inverse (n=2) over nmod %lu\n", p);
  flint_printf("---------------------------------------------------\n");
  gr_mat_t A;
  gr_mat_init(A, 2, 2, ctx);
  gr_set_ui(gr_mat_entry_ptr(A, 0, 0, ctx), 3, ctx);
  gr_set_ui(gr_mat_entry_ptr(A, 0, 1, ctx), 5, ctx);
  gr_set_ui(gr_mat_entry_ptr(A, 1, 0, ctx), 7, ctx);
  gr_set_ui(gr_mat_entry_ptr(A, 1, 1, ctx), 3, ctx);
  flint_printf("Input A:\n");
  gr_mat_print(A, ctx);
  flint_printf("\n");
  gr_mat_t G_A, H_A;
  res = gr_mat_G_H(G_A, H_A, A, ctx);
  if (res != GR_SUCCESS) {
    flint_printf("[ERROR] Could not compute generators.\n");
    gr_mat_clear(A, ctx);
    gr_ctx_clear(ctx);
    return GR_TEST_FAIL;
  }
  flint_printf("G_A: ");
  gr_mat_print(G_A, ctx);
  flint_printf("\nH_A: ");
  gr_mat_print(H_A, ctx);
  flint_printf("\n\n");
  slong rank = gr_mat_ncols(G_A, ctx);
  gr_mat_t G_D, H_D;
  gr_mat_init(G_D, 2, rank, ctx);
  gr_mat_init(H_D, 2, rank, ctx);
  res = gr_toeplitz_inverse(G_D, H_D, G_A, H_A, ctx);
  if (res != GR_SUCCESS) {
    flint_printf("[ERROR] gr_toeplitz_inverse failed.\n");
    gr_mat_clear(G_D, ctx);
    gr_mat_clear(H_D, ctx);
    gr_mat_clear(G_A, ctx);
    gr_mat_clear(H_A, ctx);
    gr_mat_clear(A, ctx);
    gr_ctx_clear(ctx);
    return GR_TEST_FAIL;
  }
  flint_printf("G_D (inv generators): ");
  gr_mat_print(G_D, ctx);
  flint_printf("\nH_D (inv generators): ");
  gr_mat_print(H_D, ctx);
  flint_printf("\n\n");

  // A^{-1} from generators
  gr_mat_t A_inv;
  gr_mat_init(A_inv, 2, 2, ctx);
  res = gr_mat_reconstruct_A_safe(A_inv, G_D, H_D, ctx);
  if (res != GR_SUCCESS) {
    flint_printf("[ERROR] Could not reconstruct A_inv.\n");
    gr_mat_clear(A_inv, ctx);
    gr_mat_clear(G_D, ctx);
    gr_mat_clear(H_D, ctx);
    gr_mat_clear(G_A, ctx);
    gr_mat_clear(H_A, ctx);
    gr_mat_clear(A, ctx);
    gr_ctx_clear(ctx);
    return GR_TEST_FAIL;
  }
  flint_printf("Reconstructed A^{-1}:\n");
  gr_mat_print(A_inv, ctx);
  flint_printf("\n\n");

  // A * A^{-1} should be identity
  gr_mat_t product;
  gr_mat_init(product, 2, 2, ctx);
  res = gr_mat_mul(product, A, A_inv, ctx);
  if (res != GR_SUCCESS) {
    flint_printf("[ERROR] Matrix multiply failed.\n");
    gr_mat_clear(product, ctx);
    gr_mat_clear(A_inv, ctx);
    gr_mat_clear(G_D, ctx);
    gr_mat_clear(H_D, ctx);
    gr_mat_clear(G_A, ctx);
    gr_mat_clear(H_A, ctx);
    gr_mat_clear(A, ctx);
    gr_ctx_clear(ctx);
    return GR_TEST_FAIL;
  }
  flint_printf("A * A^{-1}:\n");
  gr_mat_print(product, ctx);
  flint_printf("\n");
  gr_mat_t identity;
  gr_mat_init(identity, 2, 2, ctx);
  gr_mat_one(identity, ctx);
  // gr_mat_equal returns truth_t (T_TRUE / T_FALSE / T_UNKNOWN), not int
  truth_t check;
  check = gr_mat_equal(product, identity, ctx);
  if (check == T_TRUE) {
    flint_printf("[SUCCESS] A * A^{-1} = I confirmed!\n");
  } else {
    flint_printf("[FAILURE] A * A^{-1} != I\n");
    res = GR_TEST_FAIL;
  }
  gr_mat_clear(identity, ctx);
  gr_mat_clear(product, ctx);
  gr_mat_clear(A_inv, ctx);

  gr_mat_clear(G_D, ctx);
  gr_mat_clear(H_D, ctx);
  gr_mat_clear(G_A, ctx);
  gr_mat_clear(H_A, ctx);
  gr_mat_clear(A, ctx);
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
  } else if (strcmp("inverse_2x2", argv[1]) == 0) {
    ok = test_inverse_2x2();
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