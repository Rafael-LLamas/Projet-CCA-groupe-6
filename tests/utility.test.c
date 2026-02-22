#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "matrix_aux.h"
#include "random_toeplitz.h"

int test_LU_detatch() {
  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, GNMOD);
  flint_printf("\n:-------: LU Detatch on random Matrix :-------:\n");
  gr_mat_t A, LU, L, U, B, C;
  slong rank, *P;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  P = flint_malloc(sizeof(slong) * 10);
  gr_mat_init(A, 5, 10, ctx);
  gr_mat_init(C, 10, 10, ctx);
  int error = gr_mat_init_set(B, A, ctx);
  if (error != 0) {
    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);
    flint_rand_clear(state);
    flint_free(P);
    gr_ctx_clear(ctx);
    return error;
  }
  gr_mat_init(LU, 5, 10, ctx);
  gr_mat_init(L, 5, 10, ctx);
  gr_mat_init(U, 10, 10, ctx);
  error = rand_quasi_toeplitz(A, 5, 10, 0, ctx);
  if (error != 0) {
    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);
    gr_mat_clear(LU, ctx);
    gr_mat_clear(L, ctx);
    gr_mat_clear(U, ctx);
    flint_rand_clear(state);
    flint_free(P);
    gr_ctx_clear(ctx);
    return error;
  }
  flint_printf("A = \n");
  gr_mat_print(A, ctx);
  flint_printf("\n**************************************************\n");
  flint_printf("LU = \n");
  error = gr_mat_lu(&rank, P, LU, A, 0, ctx);
  if (error != 0) {
    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);
    gr_mat_clear(LU, ctx);
    gr_mat_clear(L, ctx);
    gr_mat_clear(U, ctx);
    flint_rand_clear(state);
    flint_free(P);
    gr_ctx_clear(ctx);
    return error;
  }
  gr_mat_print(LU, ctx);
  flint_printf("\n**************************************************\n");
  error = gr_mat_lu_detach(L, U, LU, ctx);
  if (error != 0) {
    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);
    gr_mat_clear(LU, ctx);
    gr_mat_clear(L, ctx);
    gr_mat_clear(U, ctx);
    flint_rand_clear(state);
    flint_free(P);
    gr_ctx_clear(ctx);
    return error;
  }
  flint_printf("U = \n");
  gr_mat_print(U, ctx);
  flint_printf("\n\n");
  flint_printf("L = \n");
  gr_mat_print(L, ctx);
  flint_printf("\n**************************************************\n");
  error = gr_mat_mul(B, L, U, ctx);
  if (error != 0) {
    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);
    gr_mat_clear(LU, ctx);
    gr_mat_clear(L, ctx);
    gr_mat_clear(U, ctx);
    flint_rand_clear(state);
    flint_free(P);
    gr_ctx_clear(ctx);
    return error;
  }
  flint_printf("B = \n");
  gr_mat_print(B, ctx);
  flint_printf("\n**************************************************\n");
  error = gr_mat_random_manual(C, state, 1009, ctx);
  if (error != 0) {
    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);
    gr_mat_clear(LU, ctx);
    gr_mat_clear(L, ctx);
    gr_mat_clear(U, ctx);
    flint_rand_clear(state);
    flint_free(P);
    gr_ctx_clear(ctx);
    return error;
  }
  flint_printf("C = \n");
  gr_mat_print(C, ctx);
  gr_mat_clear(A, ctx);
  gr_mat_clear(LU, ctx);
  gr_mat_clear(L, ctx);
  gr_mat_clear(U, ctx);
  gr_mat_clear(B, ctx);
  gr_mat_clear(C, ctx);
  flint_rand_clear(state);
  flint_free(P);
  gr_ctx_clear(ctx);
  return GR_SUCCESS;
}

int test_LU_detatch_2() {
  gr_ctx_t ctx;
  flint_printf("\n\n:-------: LU Detatch on weird Matrix :-------:\n");
  gr_mat_t A, LU, L, U, B;
  slong rank, *P;
  gr_ctx_init_nmod(ctx, GNMOD);
  P = flint_malloc(sizeof(slong) * 10);
  gr_mat_init(A, 5, 4, ctx);
  gr_mat_init(LU, 5, 4, ctx);
  gr_mat_init(L, 5, 5, ctx);
  gr_mat_init(U, 5, 4, ctx);
  gr_mat_init(B, 5, 4, ctx);
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 0, ctx), 1, ctx)); // Row 0: 1, 2, 1, 2, 1
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 1, ctx), 2, ctx));
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 2, ctx), 1, ctx));
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 3, ctx), 2, ctx));
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 4, ctx), 1, ctx));
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 0, ctx), 2, ctx)); // Row 1: 2, 4, 3, 6, 1
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 1, ctx), 4, ctx));
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 2, ctx), 3, ctx));
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 3, ctx), 6, ctx));
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 4, ctx), 1, ctx));
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 0, ctx), 3, ctx)); // Row 2: 3, 6, 3, 6, 5
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 1, ctx), 6, ctx));
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 2, ctx), 3, ctx));
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 3, ctx), 6, ctx));
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 4, ctx), 5, ctx));
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 3, 0, ctx), 4, ctx)); // Row 3: 4, 8, 4, 8, 6
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 3, 1, ctx), 8, ctx));
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 3, 2, ctx), 4, ctx));
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 3, 3, ctx), 8, ctx));
  FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 3, 4, ctx), 6, ctx));
  // FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 4, 0, ctx), 5, ctx)); // Row 4:  5, 10, 5, 10, 5,
  // FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 4, 1, ctx), 10, ctx));
  // FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 4, 2, ctx), 5, ctx));
  // FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 4, 3, ctx), 10, ctx));
  // FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 4, 4, ctx), 5, ctx));
  flint_printf("A = \n");
  gr_mat_print(A, ctx);
  flint_printf("\n**************************************************\n");
  flint_printf("LU = \n");
  FLINT_CHECK(gr_mat_lu(&rank, P, LU, A, 0, ctx));
  gr_mat_print(LU, ctx);
  flint_printf("\n**************************************************\n");
  FLINT_CHECK(gr_mat_lu_detach(L, U, LU, ctx));
  flint_printf("U = \n");
  gr_mat_print(U, ctx);
  flint_printf("\n\n");
  flint_printf("L = \n");
  gr_mat_print(L, ctx);
  flint_printf("\n**************************************************\n");
  FLINT_CHECK(gr_mat_mul(B, L, U, ctx));
  flint_printf("B = \n");
  gr_mat_print(B, ctx);
  flint_printf("\n");
  // for (slong i = 0; i < gr_mat_nrows(B, ctx); i++) FLINT_CHECK(gr_mat_move_row(B, i, P[i], ctx)); // this fails
  flint_printf("B but (not yet) permutated = \n");
  gr_mat_print(B, ctx);
  flint_printf("\n");
  gr_mat_clear(A, ctx);
  gr_mat_clear(LU, ctx);
  gr_mat_clear(L, ctx);
  gr_mat_clear(U, ctx);
  gr_mat_clear(B, ctx);
  flint_free(P);
  gr_ctx_clear(ctx);
  return GR_SUCCESS;
}

int main(int argc, char *argv[]) {
  if (argc == 1) { return GR_UNABLE; }
  // start test
  fprintf(stderr, "=> Start test \"%s\"\n", argv[1]);
  int ok = GR_SUCCESS;
  if (strcmp("LU_detatch", argv[1]) == 0) {
    ok = test_LU_detatch();
  } else if (strcmp("LU_detatch2", argv[1]) == 0) {
    ok = test_LU_detatch_2();
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