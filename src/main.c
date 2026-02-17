#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "displacement_matrices.h"
#include "operations/addition.h"
#include "operations/multiplication.h"
#include "random_toeplitz.h"
#include "utility/matrix_aux.h"

/*

Example piece of code to execute.
The tests are moved to project/tests,
$ make test
to test out. Check out LastTest.log
in Testing/Temporary

*/

int main() {
  gr_ctx_t ctx;
  gr_ctx_init_nmod(ctx, GNMOD);
  flint_rand_t state;
  flint_rand_init(state);
  slong N = 5;
  gr_mat_t A, B, C_ref, C_gen; // Matrix Declarations
  gr_mat_init(A, N, N, ctx);
  gr_mat_init(B, N, N, ctx);
  gr_mat_init(C_ref, N, N, ctx);
  gr_mat_init(C_gen, N, N, ctx);
  gr_mat_t GA, HA, GB, HB, GC, HC; // Generator Declarations

  flint_printf("--- 1. Generating Random Quasi-Toeplitz Matrices ---\n");
  rand_quasi_toeplitz(A, N, N, 2, ctx);
  rand_quasi_toeplitz(B, N, N, 4, ctx);
  flint_printf("Matrix A:\n");
  gr_mat_print(A, ctx);
  flint_printf("\nMatrix B:\n");
  gr_mat_print(B, ctx);
  flint_printf("\n");

  flint_printf("--- 2. Computing Generators (Displacement Decomposition) ---\n");
  gr_mat_G_H(GA, HA, A, ctx);
  gr_mat_G_H(GB, HB, B, ctx);
  slong rank_A = gr_mat_ncols(GA, ctx);
  slong rank_B = gr_mat_ncols(GB, ctx);
  flint_printf("Displacement Rank of A: %ld\n", rank_A);
  flint_printf("Displacement Rank of B: %ld\n", rank_B);
  gr_mat_t A_rec;
  gr_mat_init(A_rec, N, N, ctx);
  gr_mat_reconstruct_A_safe(A_rec, GA, HA, ctx);
  if (gr_mat_equal(A, A_rec, ctx) == T_TRUE) {
    flint_printf("[Check] Reconstruction of A from generators successful.\n");
  } else {
    flint_printf("[Check] Reconstruction of A FAILED.\n");
  }
  gr_mat_clear(A_rec, ctx);
  flint_printf("\n");

  flint_printf("--- 3. Addition via Generators ---\n");
  FLINT_CHECK(gr_mat_add(C_ref, A, B, ctx));
  gr_mat_addition_generateur(GA, HA, GB, HB, GC, HC, ctx);
  FLINT_CHECK(gr_mat_zero(C_gen, ctx));
  gr_mat_reconstruct_A_safe(C_gen, GC, HC, ctx);
  flint_printf("Rank of Sum (A+B): %ld\n", gr_mat_ncols(GC, ctx));
  if (gr_mat_equal(C_ref, C_gen, ctx) == T_TRUE) {
    flint_printf("[Success] Generator Addition matches Standard Addition.\n");
  } else {
    flint_printf("[Failure] Generator Addition mismatch.\n");
  }
  flint_printf("\n");

  flint_printf("--- 4. Multiplication via Generators TODO ---\n");

  gr_mat_clear(A, ctx);
  gr_mat_clear(B, ctx);
  gr_mat_clear(C_ref, ctx);
  gr_mat_clear(C_gen, ctx);
  gr_mat_clear(GA, ctx);
  gr_mat_clear(HA, ctx);
  gr_mat_clear(GB, ctx);
  gr_mat_clear(HB, ctx);
  gr_mat_clear(GC, ctx);
  gr_mat_clear(HC, ctx);
  flint_rand_clear(state);
  gr_ctx_clear(ctx);
  return GR_SUCCESS;
}