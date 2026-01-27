#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include <stdlib.h>
#include <time.h>

#define FLINT_CHECK(x)                                                                                                 \
  do {                                                                                                                 \
    int _status = (x);                                                                                                 \
    if (_status != GR_SUCCESS) return _status;                                                                         \
  } while (0)

/*
We are currently focusing on square Toeplitz matrices.
We can represent a Toeplitz matrix with only 2 vectors,
thus reducing its size + time of calculation.

This file is the first step in calculation of this
displacement matrix (∇A):

∇A = A − ZAZ^T

This essentially shifts the matrix A one step to the right
and bottom, leaving the first line and the first column
with 0s, minus A, creates 2 vectors representing the
Toeplitz matrix.

Some vocabulary:
    Z   : Lower Shift
    Z^T : Upper Shift

Later we will focus on non square toeplitz matrices,
block toeplitz matrices.

Structure for optimized version:
We know that the displacemenet matrix is an L shaped one.

                    L U U U ... U
                    L 0 0 0 ... 0
                    L 0 0 0 ... 0
                    ... 0 0 ... 0
                    L 0 0 0 ... 0

Which is 2 vectors. For future proofing, we will store them
in 2 matrices.

        L   1
        L   0               1 0 0 0 ... 0
        L   0       and     0 U U U ... U
        ... 0
        L   0

We will name these two matrices G and H^T as marked in
Definition 1.6
For different rangs, the size of these will change.
Doing this way allows us to use gr_mat_mul directly rather
than translating two instances of vectors.
*/

/**
 * @brief This is the safe-book version relying purely on the optimisation of FLINT calculations.
 * Complexity: O(n^3) with naive, but this is flint so probably less.
 * @param[out] D Resulting L shaped matrix of xnx
 * @param[in] A Input square nxn toeplitz matrix
 */
int gr_mat_displacement_square(gr_mat_t D, gr_mat_t A, gr_ctx_t ctx) {
  slong n = gr_mat_nrows(A, ctx);
  gr_mat_t Z, ZT, Temp1, Temp2;
  gr_mat_init(Z, n, n, ctx);
  gr_mat_init(ZT, n, n, ctx);
  gr_mat_init(Temp1, n, n, ctx);
  gr_mat_init(Temp2, n, n, ctx);
  FLINT_CHECK(gr_mat_zero(Z, ctx));
  gr_ptr val;
  for (slong i = 1; i < n; i++) { // create Z
    val = gr_mat_entry_ptr(Z, i, i - 1, ctx);
    FLINT_CHECK(gr_one(val, ctx));
  }
  FLINT_CHECK(gr_mat_transpose(ZT, Z, ctx));     // Z^T
  FLINT_CHECK(gr_mat_mul(Temp1, A, ZT, ctx));    // Temp1 = A * Z^T
  FLINT_CHECK(gr_mat_mul(Temp2, Z, Temp1, ctx)); // Temp2 = Z * Temp1  (which is Z * A * Z^T)
  FLINT_CHECK(gr_mat_sub(D, A, Temp2, ctx));     // D = A - Temp2
  gr_mat_clear(Z, ctx);
  gr_mat_clear(ZT, ctx);
  gr_mat_clear(Temp1, ctx);
  gr_mat_clear(Temp2, ctx);
  return GR_SUCCESS;
}

/**
 * @brief Computes the displacement generators G and H for a nxn Toeplitz matrix A.
 *  * Complexity: O(2n) ->
 * Instead of performing in full arithmetic, we copy the first column
 * and first row of A.
 * @param[out] G Output matrix of size nx2.
 * Must be initialized before calling.
 * Will contain the vertical vectors.
 * @param[out] H Output matrix of size nx2.
 * Must be initialized before calling.
 * Will contain the horizontal generator vectors (stored as columns).
 * @param[in] A The input square Toeplitz matrix of size nxn.
 * @param[in] ctx The FLINT generic ring context.
 * @return GR_SUCCESS on success, or GR_DOMAIN if dimensions do not match.
 */
int gr_mat_displacement_square_opt(gr_mat_t G, gr_mat_t H, gr_mat_t A, gr_ctx_t ctx) {
  slong n = gr_mat_nrows(A, ctx); // nxn of A.
  if (gr_mat_nrows(G, ctx) != n || gr_mat_ncols(G, ctx) != 2 || gr_mat_nrows(H, ctx) != n || gr_mat_ncols(H, ctx) != 2)
    return GR_DOMAIN;               // safety check
  FLINT_CHECK(gr_mat_zero(G, ctx)); // set G and H to 0 for safety
  FLINT_CHECK(gr_mat_zero(H, ctx));

  gr_ptr src_val, dst_val; // pointers

  for (slong i = 0; i < n; i++) { // G gets the first column of A
    src_val = gr_mat_entry_ptr(A, i, 0, ctx);
    dst_val = gr_mat_entry_ptr(G, i, 0, ctx);
    FLINT_CHECK(gr_set(dst_val, src_val, ctx));
  }
  dst_val = gr_mat_entry_ptr(G, 0, 1, ctx);
  FLINT_CHECK(gr_one(dst_val, ctx)); // G[0][1] = 1

  for (slong j = 1; j < n; j++) { // H gets the first row of A
    src_val = gr_mat_entry_ptr(A, 0, j, ctx);
    dst_val = gr_mat_entry_ptr(H, j, 1, ctx);
    FLINT_CHECK(gr_set(dst_val, src_val, ctx));
  }
  dst_val = gr_mat_entry_ptr(H, 0, 0, ctx);
  FLINT_CHECK(gr_one(dst_val, ctx)); // H[0][0] = 1

  return GR_SUCCESS;
}

// temporary testing
int main() {
  gr_ctx_t ctx;
  gr_ctx_init_fmpz(ctx);
  slong n = 500;
  flint_printf("Matrix Size: %wd x %wd\n", n, n);
  gr_mat_t A, G, H, HT, Product, D_slow;
  gr_mat_init(A, n, n, ctx);
  gr_mat_init(G, n, 2, ctx);
  gr_mat_init(H, n, 2, ctx);
  gr_mat_init(HT, 2, n, ctx);
  gr_mat_init(Product, n, n, ctx);
  gr_mat_init(D_slow, n, n, ctx);
  // this is a temporary way to create a toeplitz matrix dont mind this its not random
  void *val = gr_heap_init(ctx);
  gr_ptr entry_ptr;
  for (slong i = 0; i < n; i++) {
    for (slong j = 0; j < n; j++) {
      (void)gr_set_si(val, i - j, ctx);
      entry_ptr = gr_mat_entry_ptr(A, i, j, ctx);
      (void)gr_set(entry_ptr, val, ctx);
    }
  }
  gr_heap_clear(val, ctx);

  // Test 1: Safe Book
  clock_t start = clock();
  if (gr_mat_displacement_square(D_slow, A, ctx) != GR_SUCCESS) {
    flint_printf("Safe book method failed.\n");
    return EXIT_FAILURE;
  }
  clock_t end = clock();
  double time_slow = (double)(end - start) / CLOCKS_PER_SEC;
  flint_printf("\nArithmetic Method: %.6f seconds\n", time_slow);

  // Test 2: Copying
  start = clock();
  if (gr_mat_displacement_square_opt(G, H, A, ctx) != GR_SUCCESS) {
    flint_printf("Copying method failed.\n");
    return EXIT_FAILURE;
  }
  end = clock();
  double time_fast = (double)(end - start) / CLOCKS_PER_SEC;
  flint_printf("Copying Method: %.6f seconds\n", time_fast);

  FLINT_CHECK(gr_mat_transpose(HT, H, ctx));
  FLINT_CHECK(gr_mat_mul(Product, G, HT, ctx));

  truth_t is_equal = gr_mat_equal(Product, D_slow, ctx);
  flint_printf("\n------------------ Result: --------------------\n");
  if (is_equal == T_TRUE)
    flint_printf("[SUCCESS] The Copying Method result matches the Arithmetic Method exactly so they must be correct!!! (unless.....)\n\n");
  else
    flint_printf("[FAILURE] The results are different. Something is wrong D:\n\n");

  gr_mat_clear(A, ctx);
  gr_mat_clear(G, ctx);
  gr_mat_clear(H, ctx);
  gr_mat_clear(HT, ctx);
  gr_mat_clear(Product, ctx);
  gr_mat_clear(D_slow, ctx);
  gr_ctx_clear(ctx);
  return EXIT_SUCCESS;
}