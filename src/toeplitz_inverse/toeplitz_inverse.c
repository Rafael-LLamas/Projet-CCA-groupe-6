#include "addition.h"
#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "multiplication.h"

/*
Algorithm 10.1 Algorithme de type Strassen pour inverser une matrice
quasi-Toeplitz, from [1].

Only square and non singular matrices are invertible. We will reject
any input of a non square matrix.
*/

int gr_toeplitz_inverse(gr_mat_t G_D, gr_mat_t H_D, gr_mat_t G_A, gr_mat_t H_A, gr_ctx_t ctx) {
  int res = GR_SUCCESS;
  // TODO: can maybe add a check of its determinant once implemented the toeplitz version (det must not be zero)
  gr_mat_zero(G_D, ctx);
  gr_mat_zero(H_D, ctx);
  slong rank = gr_mat_ncols(G_A, ctx);

  // 1: base case if n = 1, return A^{-1}
  if (gr_mat_nrows(G_A, ctx) == 1) {
    // sum(Ga[0, 0:rank] * Ha[0, 0:rank])
    for (slong i = 0; i < rank; i++) {
      // Ga[0,i] * Ha[0,i] -(temp)-> H_D[0,0], then -> G_D[0, 0]
      gr_mul(gr_mat_entry_ptr(H_D, 0, 0, ctx), gr_mat_entry_ptr(G_A, 0, i, ctx), gr_mat_entry_ptr(H_A, 0, i, ctx), ctx);
      gr_add(gr_mat_entry_ptr(G_D, 0, 0, ctx), gr_mat_entry_ptr(G_D, 0, 0, ctx), gr_mat_entry_ptr(H_D, 0, 0, ctx), ctx);
    }

    // take the inv value in nmod (stored in G_D[0,0])
    res = gr_inv(gr_mat_entry_ptr(G_D, 0, 0, ctx), gr_mat_entry_ptr(G_D, 0, 0, ctx), ctx);
    if (res != GR_SUCCESS) { return res; }

    gr_one(gr_mat_entry_ptr(H_D, 0, 0, ctx), ctx);
    return res;
  }

  // 2: calculate the generators a b c d for n/2
  // a: Gtop Htop | b: Gtop Hbottom
  // c: Gbottom Htop | d: Gbottom Hbottom
  gr_mat_t G_top, H_top, G_bottom, H_bottom;
  slong n1 = (gr_mat_nrows(G_A, ctx) + 1) / 2;
  slong n2 = gr_mat_nrows(G_A, ctx) / 2;
  gr_mat_window_init(G_top, G_A, 0, 0, n1, rank, ctx);                         // Top half of G
  gr_mat_window_init(H_top, H_A, 0, 0, n1, rank, ctx);                         // Top half of H
  gr_mat_window_init(G_bottom, G_A, n1, 0, rank, gr_mat_ncols(G_A, ctx), ctx); // bottom half of G
  gr_mat_window_init(H_bottom, H_A, n1, 0, rank, gr_mat_ncols(G_A, ctx), ctx); // bottom half of H

  // 3: calculate rec e := a^{-1}
  gr_mat_t G_e, H_e;
  gr_mat_init(G_e, n1, gr_mat_ncols(G_A, ctx), ctx);
  gr_mat_init(H_e, n1, gr_mat_ncols(G_A, ctx), ctx);
  res = gr_toeplitz_inverse(G_e, H_e, G_top, H_top, ctx);
  if (res != GR_SUCCESS) {
    gr_mat_clear(G_e, ctx);
    gr_mat_clear(H_e, ctx);
    gr_mat_window_clear(G_top, ctx);
    gr_mat_window_clear(H_top, ctx);
    gr_mat_window_clear(G_bottom, ctx);
    gr_mat_window_clear(H_bottom, ctx);
    return res;
  }

  // 4: calculate generators for S := d - ceb

  // 5: calculate inv generators for t := S^{-1}

  // 6: return the inv generators for A^{-1} x y z t with the formules of strassen y:= -ebt, z := -tce, x := e + ebtce

  gr_mat_window_clear(G_top, ctx);
  gr_mat_window_clear(H_top, ctx);
  gr_mat_window_clear(G_bottom, ctx);
  gr_mat_window_clear(H_bottom, ctx);
  gr_mat_clear(G_e, ctx);
  gr_mat_clear(H_e, ctx);
  return res;
}