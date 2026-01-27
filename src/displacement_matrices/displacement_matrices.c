#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include <stdlib.h>

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

Structure:
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
int gr_mat_displacement_square(gr_mat_t G, gr_mat_t H, gr_mat_t A, gr_ctx_t ctx) {
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
  gr_ctx_init_fmpz(ctx); // Integers
  slong n = 5;

  gr_mat_t A, G, H, HT, Product;
  gr_mat_init(A, n, n, ctx);
  gr_mat_init(G, n, 2, ctx);
  gr_mat_init(H, n, 2, ctx);
  gr_mat_init(HT, 2, n, ctx);
  gr_mat_init(Product, n, n, ctx);

  // temporary way to create a toeplitz dont mind this
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
  flint_printf("Matrix A:\n");
  gr_mat_print(A, ctx);

  if (gr_mat_displacement_square(G, H, A, ctx) == GR_SUCCESS) {
    flint_printf("\n Result G:\n");
    gr_mat_print(G, ctx);
    flint_printf("\n Result H:\n");
    gr_mat_print(H, ctx);
    FLINT_CHECK(gr_mat_transpose(HT, H, ctx));
    FLINT_CHECK(gr_mat_mul(Product, G, HT, ctx));
    flint_printf("\n Product (G * H^T):\n");
    gr_mat_print(Product, ctx);
    flint_printf("\n");
  } else {
    flint_printf("Error: Calculation failed.\n");
  }

  // Cleanup
  gr_mat_clear(A, ctx);
  gr_mat_clear(G, ctx);
  gr_mat_clear(H, ctx);
  gr_mat_clear(HT, ctx);
  gr_mat_clear(Product, ctx);
  gr_ctx_clear(ctx);

  return EXIT_SUCCESS;
}