#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"

#include <stdlib.h>
#include <time.h>

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

Structure for optimized storage:
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

/*
Keeping this function in mind,
Let's turn our focus to a more general case.

For any matrix A of size nxm, we find it's ∇A.
(rows x cols)

I noticed, doing on paper, its just a substraction
of diagonal values. We can optimize this operation while
evading temporary (any Z) matrices.
Instead of O(n^3), we can just iterate through n * m
elements (first n-1 * m-1 then n + m), giving us a
complexity of O(nm). Of matrix A:

                    1 3 5
                    8 7 8
                    9 8 7
                    2 1 1

We can substract the top left diagonal values starting from
(n,m) to (1,1). Resulting in a correct shift + substract
operation, giving us ∇A:

                    1 3 5
                    8 6 5
                    9 0 0
                    2 -8 -7
*/

int gr_mat_displacement(gr_mat_t D, gr_mat_t A, gr_ctx_t ctx) {
  slong n = gr_mat_nrows(A, ctx);
  slong m = gr_mat_ncols(A, ctx);
  if (gr_mat_nrows(D, ctx) != n || gr_mat_ncols(D, ctx) != m) return GR_UNABLE;

  gr_ptr ptr_cur, ptr_prev, ptr_dest;
  // inner part
  for (slong i = n - 1; i > 0; i--) {
    for (slong j = m - 1; j > 0; j--) {
      ptr_cur = gr_mat_entry_ptr(A, i, j, ctx);
      ptr_prev = gr_mat_entry_ptr(A, i - 1, j - 1, ctx);
      ptr_dest = gr_mat_entry_ptr(D, i, j, ctx);
      FLINT_CHECK(gr_sub(ptr_dest, ptr_cur, ptr_prev, ctx)); // D[i,j] = A[i,j] - A[i-1, j-1]
    }
  }
  // the rest
  for (slong i = 0; i < n; i++)
    FLINT_CHECK(gr_set(gr_mat_entry_ptr(D, i, 0, ctx), gr_mat_entry_ptr(A, i, 0, ctx), ctx));
  for (slong j = 1; j < m; j++)
    FLINT_CHECK(gr_set(gr_mat_entry_ptr(D, 0, j, ctx), gr_mat_entry_ptr(A, 0, j, ctx), ctx));
  return GR_SUCCESS;
}

/*
Now that we can corectly distinguish where we can optimize.

Any matrix call from gr_mat_displacement is going to return ∇A,
for the previous example, row 2 (9 0 0) we have 2 zeros. It seems
small for this matrix in particular but if it was a big matrix,
we would have gone a lot of 0s.

We can see this as the ones marked form a 2x2 Toeplitz matrix:

                    1 3 5
                    8 *7 *8
                    9 *8 *7
                    2 1 1


*/

int test_displacement_matrices() {
  flint_printf("*----------* Displacement Matrix Test *----------*\n");
  gr_ctx_t ctx;
  gr_ctx_init_fmpz(ctx);

  flint_printf(":-------: Manual nxn Toeplitz Matrix Deplacement Test :-------:\n");
  {
    gr_mat_t A, D;
    gr_mat_init(A, 3, 3, ctx);
    gr_mat_init(D, 3, 3, ctx);
    // feast your eyes
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 0, ctx), 1, ctx)); // Row 0: 1, 3, 4
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 1, ctx), 3, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 2, ctx), 4, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 0, ctx), 8, ctx)); // Row 1: 8, 1, 3
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 1, ctx), 1, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 2, ctx), 3, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 0, ctx), 9, ctx)); // Row 2: 9, 8, 1
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 1, ctx), 8, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 2, ctx), 1, ctx));
    flint_printf("Input Manual Matrix A:\n");
    gr_mat_print(A, ctx);
    flint_printf("\n");
    gr_mat_displacement(D, A, ctx);
    flint_printf("Displacement Matrix (∇A):\n");
    gr_mat_print(D, ctx);
    flint_printf("\n");
    gr_mat_clear(A, ctx);
    gr_mat_clear(D, ctx);
  }

  flint_printf(":-------: Manual nxm Matrix Test :-------:\n");
  {
    gr_mat_t A, D;
    gr_mat_init(A, 4, 3, ctx);
    gr_mat_init(D, 4, 3, ctx);
    // (i am too tired to do a double triangular loop future me if ur reading this... dew it -palpatine)
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 0, ctx), 1, ctx)); // Row 0: 1, 3, 4
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 1, ctx), 3, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 0, 2, ctx), 4, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 0, ctx), 8, ctx)); // Row 1: 8, 7, 8
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 1, ctx), 7, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 1, 2, ctx), 8, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 0, ctx), 9, ctx)); // Row 2: 9, 8, 7
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 1, ctx), 8, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 2, 2, ctx), 7, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 3, 0, ctx), 2, ctx)); // Row 3: 2, 1, 1
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 3, 1, ctx), 1, ctx));
    FLINT_CHECK(gr_set_si(gr_mat_entry_ptr(A, 3, 2, ctx), 1, ctx));
    flint_printf("Input Manual Matrix A:\n");
    gr_mat_print(A, ctx);
    flint_printf("\n");
    gr_mat_displacement(D, A, ctx);
    flint_printf("Displacement Matrix (∇A):\n");
    gr_mat_print(D, ctx);
    flint_printf("\n");
    slong rank;
    slong *P = flint_malloc(4 * sizeof(slong));
    gr_mat_t LU;
    gr_mat_init(LU, 4, 3, ctx);
    if (gr_mat_lu(&rank, P, LU, D, 0, ctx) == GR_SUCCESS) {
      flint_printf("LU Decomposition successful.\n");
      gr_mat_print(LU, ctx);
    } else
      flint_printf("LU Decomposition failed.\n"); // due to maths, 1.5 is not an integer
    flint_free(P);
    gr_mat_clear(LU, ctx);
    gr_mat_clear(A, ctx);
    gr_mat_clear(D, ctx);
  }
  gr_ctx_clear(ctx);
  return GR_SUCCESS;
}