#include "addition.h"
#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "multiplication.h"

/*
Algorithm 10.1 Algorithme de type Strassen pour inverser une matrice
quasi-Toeplitz from [1].

Only square and non singular matrices are invertible. We will reject
any input of a non square matrix.
*/

int gr_toeplitz_inverse(gr_mat_t G_D, gr_mat_t H_D, gr_mat_t G_A, gr_mat_t H_A, gr_ctx_t ctx) {
  int res = GR_SUCCESS;
  // TODO: add a check of its determinant once implemented the toeplitz version (det must not be zero)

  // base case
  if (gr_mat_nrows(G_A, ctx) == 1) {
    gr_zero(gr_mat_entry_ptr(G_D, 0, 0, ctx), ctx);

    // sum(Ga[0, 0:rank] * Ha[0, 0:rank])
    for (slong i = 0; i < gr_mat_ncols(G_A, ctx); i++) {
      // Ga[0,i] * Ha[0,i] -(temp)-> H_D[0,0], then -> G_D[0, 0]
      gr_mul(gr_mat_entry_ptr(H_D, 0, 0, ctx), gr_mat_entry_ptr(G_A, 0, i, ctx), gr_mat_entry_ptr(H_A, 0, i, ctx), ctx);
      gr_add(gr_mat_entry_ptr(G_D, 0, 0, ctx), gr_mat_entry_ptr(G_D, 0, 0, ctx), gr_mat_entry_ptr(H_D, 0, 0, ctx), ctx);
    }

    // take the inv value in nmod (stored in G_D[0,0])
    res = gr_inv(gr_mat_entry_ptr(G_D, 0, 0, ctx), gr_mat_entry_ptr(G_D, 0, 0, ctx), ctx);
    if (res != GR_SUCCESS) { return res; }

    gr_one(gr_mat_entry_ptr(H_D, 0, 0, ctx), ctx);

    // if rank > 2 better to put everything else to 0
    for (slong i = 2; i < gr_mat_ncols(G_A, ctx); i++) {
      gr_zero(gr_mat_entry_ptr(G_D, 0, i, ctx), ctx);
      gr_zero(gr_mat_entry_ptr(H_D, 0, i, ctx), ctx);
    }
    return res;
  }

  // a b c d

  //   gr_mat_t G_a, H_a;
  //   gr_mat_window_init(G_a, G_A, 0, 0, n / 2, rank, ctx); // Top half of G
  //   gr_mat_window_init(H_a, H_A, 0, 0, n / 2, rank, ctx); // Top half of H

  //   gr_toeplitz_inverse_strassen(G_e, H_e, G_a, H_a, ctx);

  //   gr_mat_window_clear(G_a, ctx);
  //   gr_mat_window_clear(H_a, ctx);
  return res;
}

// for debug
int main() {
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
  int res = gr_toeplitz_inverse(G_D, H_D, G_A, H_A, ctx);
  if (res == GR_SUCCESS) {
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
  }
  gr_mat_clear(G_A, ctx);
  gr_mat_clear(H_A, ctx);
  gr_mat_clear(G_D, ctx);
  gr_mat_clear(H_D, ctx);
  gr_ctx_clear(ctx);
  return 0;
}