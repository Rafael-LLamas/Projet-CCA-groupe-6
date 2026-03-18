#include "addition.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "multiplication.h"

/*
Algorithm 10.1 Algorithme de type Strassen pour inverser une matrice
quasi-Toeplitz, from [1].
*/

int gr_toeplitz_inverse(gr_mat_t G_D, gr_mat_t H_D, gr_mat_t G_A, gr_mat_t H_A, gr_ctx_t ctx) {

  /*
   * The Strassed inversion algorithm adapted to the generators of Toeplitz Matrices
   *
   *
   */

  int status = GR_SUCCESS;
  slong n = gr_mat_nrows(G_A, ctx);
  slong rank = gr_mat_ncols(G_A, ctx);
  // TODO: can maybe add a check of its determinant once implemented the toeplitz version (det must not be zero)
  status |= gr_mat_zero(G_D, ctx);
  status |= gr_mat_zero(H_D, ctx);

  // base case n = 1
  if (n == 1) {
    if (rank == 0) { return GR_UNABLE; }
    status |= gr_inv(gr_mat_entry_ptr(G_D, 0, 0, ctx), gr_mat_entry_ptr(G_A, 0, 0, ctx), ctx);
    status |= gr_inv(gr_mat_entry_ptr(H_D, 0, 0, ctx), gr_mat_entry_ptr(H_A, 0, 0, ctx), ctx);
    return status;
  }

  // 2: calculate the generators a b c d for n/2
  // a: Gtop Htop | b: Gtop Hbottom
  // c: Gbottom Htop | d: (Gbottom Hbottom) + a
  slong n1 = (n + 1) / 2;
  slong n2 = n / 2;

  gr_mat_t G_top, H_top, G_bottom, H_bottom;
  gr_mat_window_init(G_top, G_A, 0, 0, n1, rank, ctx);
  gr_mat_window_init(H_top, H_A, 0, 0, n1, rank, ctx);
  gr_mat_window_init(G_bottom, G_A, n1, 0, n, rank, ctx);
  gr_mat_window_init(H_bottom, H_A, n1, 0, n, rank, ctx);

  // 3: calculate rec e := a^{-1}
  gr_mat_t G_e, H_e;
  gr_mat_init(G_e, n1, rank, ctx);
  gr_mat_init(H_e, n1, rank, ctx);
  status |= gr_toeplitz_inverse(G_e, H_e, G_top, H_top, ctx);
  if (status != GR_SUCCESS) goto free_e_top_bottom;

  // 4: calculate generators for S := d - ceb

  // 4.1 ce = c * e
  gr_mat_t G_ce, H_ce;
  status |= gr_mat_mul_generator(G_ce, H_ce, G_bottom, H_top, G_e, H_e, ctx);
  if (status != GR_SUCCESS) {
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_e_top_bottom;
  }

  // 4.2 eb = e * b
  gr_mat_t G_eb, H_eb;
  status |= gr_mat_mul_generator(G_eb, H_eb, G_e, H_e, G_top, H_bottom, ctx);
  if (status != GR_SUCCESS) {
    gr_mat_clear(G_eb, ctx);
    gr_mat_clear(H_eb, ctx);
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_e_top_bottom;
  }

  // 4.3 ceb = ce * b, then negate G in-place
  gr_mat_t G_ceb, H_ceb;
  status |= gr_mat_mul_generator(G_ceb, H_ceb, G_ce, H_ce, G_top, H_bottom, ctx);
  status |= gr_mat_neg(G_ceb, G_ceb, ctx);
  if (status != GR_SUCCESS) {
    gr_mat_clear(G_ceb, ctx);
    gr_mat_clear(H_ceb, ctx);
    gr_mat_clear(G_eb, ctx);
    gr_mat_clear(H_eb, ctx);
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_e_top_bottom;
  }

  // 4.4 S = d + (−ceb)
  gr_mat_t G_S, H_S;
  status |= gr_mat_addition_generateur(G_bottom, H_bottom, G_ceb, H_ceb, G_S, H_S, ctx);
  gr_mat_clear(G_ceb, ctx);
  gr_mat_clear(H_ceb, ctx);
  if (status != GR_SUCCESS) {
    gr_mat_clear(G_S, ctx);
    gr_mat_clear(H_S, ctx);
    gr_mat_clear(G_eb, ctx);
    gr_mat_clear(H_eb, ctx);
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_e_top_bottom;
  }

  // 5: calculate inv generators for t := S^{-1}
  gr_mat_t G_t, H_t;
  gr_mat_init(G_t, n2, gr_mat_ncols(G_S, ctx), ctx);
  gr_mat_init(H_t, n2, gr_mat_ncols(H_S, ctx), ctx);
  status |= gr_toeplitz_inverse(G_t, H_t, G_S, H_S, ctx);
  gr_mat_clear(G_S, ctx);
  gr_mat_clear(H_S, ctx);
  if (status != GR_SUCCESS) {
    gr_mat_clear(G_eb, ctx);
    gr_mat_clear(H_eb, ctx);
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_t;
  }

  // 6: return the inv generators for A^{-1} x y z t with the formules of strassen y:= -ebt, z := -tce, x := e + ebtce

  // ebt = eb * t
  gr_mat_t G_ebt, H_ebt;
  status |= gr_mat_mul_generator(G_ebt, H_ebt, G_eb, H_eb, G_t, H_t, ctx);
  gr_mat_clear(G_eb, ctx);
  gr_mat_clear(H_eb, ctx);
  if (status != GR_SUCCESS) {
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_ebt;
  }

  // tce = t * ce
  gr_mat_t G_tce, H_tce;
  status |= gr_mat_mul_generator(G_tce, H_tce, G_t, H_t, G_ce, H_ce, ctx);
  if (status != GR_SUCCESS) {
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_tce;
  }

  // ebtce = ebt * ce
  gr_mat_t G_ebtce, H_ebtce;
  status |= gr_mat_mul_generator(G_ebtce, H_ebtce, G_ebt, H_ebt, G_ce, H_ce, ctx);
  gr_mat_clear(G_ce, ctx);
  gr_mat_clear(H_ce, ctx);
  if (status != GR_SUCCESS) {
    gr_mat_clear(G_ebtce, ctx);
    gr_mat_clear(H_ebtce, ctx);
    goto free_tce;
  }

  // x = e + ebtce
  gr_mat_t G_x, H_x;
  status |= gr_mat_addition_generateur(G_e, H_e, G_ebtce, H_ebtce, G_x, H_x, ctx);
  gr_mat_clear(G_ebtce, ctx);
  gr_mat_clear(H_ebtce, ctx);
  if (status != GR_SUCCESS) { goto free_x; }

  // y = -ebt (in-place negate G_ebt)
  status |= gr_mat_neg(G_ebt, G_ebt, ctx);
  if (status != GR_SUCCESS) { goto free_x; }

  // z = -tce (in-place negate G_tce)
  status |= gr_mat_neg(G_tce, G_tce, ctx);
  if (status != GR_SUCCESS) { goto free_x; }

  // G_D / H_D
  slong rx = gr_mat_ncols(G_x, ctx);
  slong ry = gr_mat_ncols(G_ebt, ctx);
  slong rz = gr_mat_ncols(G_tce, ctx);
  slong rt = gr_mat_ncols(G_t, ctx);
  slong total_rank = rx + ry + rz + rt;

  gr_mat_clear(G_D, ctx);
  gr_mat_clear(H_D, ctx);
  gr_mat_init(G_D, n, total_rank, ctx);
  gr_mat_init(H_D, n, total_rank, ctx);
  status |= gr_mat_zero(G_D, ctx);
  status |= gr_mat_zero(H_D, ctx);

  const gr_mat_struct *G_blk[4] = {G_x, G_ebt, G_tce, G_t};
  const gr_mat_struct *H_blk[4] = {H_x, H_ebt, H_tce, H_t};
  const slong G_rs[4] = {0, 0, n1, n1};
  const slong H_rs[4] = {0, n1, 0, n1};
  const slong blk_rows[4] = {n1, n1, n2, n2};
  const slong col_off[4] = {0, rx, rx + ry, rx + ry + rz};

  for (slong b = 0; b < 4; b++) {
    slong rk = gr_mat_ncols(G_blk[b], ctx);
    for (slong k = 0; k < rk; k++) {
      slong col = col_off[b] + k;
      for (slong i = 0; i < blk_rows[b]; i++) {
        status |= gr_set(gr_mat_entry_ptr(G_D, G_rs[b] + i, col, ctx), gr_mat_entry_srcptr(G_blk[b], i, k, ctx), ctx);
        status |= gr_set(gr_mat_entry_ptr(H_D, H_rs[b] + i, col, ctx), gr_mat_entry_srcptr(H_blk[b], i, k, ctx), ctx);
      }
    }
  }

free_x:
  gr_mat_clear(G_x, ctx);
  gr_mat_clear(H_x, ctx);
free_tce:
  gr_mat_clear(G_tce, ctx);
  gr_mat_clear(H_tce, ctx);
free_ebt:
  gr_mat_clear(G_ebt, ctx);
  gr_mat_clear(H_ebt, ctx);
free_t:
  gr_mat_clear(G_t, ctx);
  gr_mat_clear(H_t, ctx);
free_e_top_bottom:
  gr_mat_clear(G_e, ctx);
  gr_mat_clear(H_e, ctx);
  gr_mat_window_clear(G_top, ctx);
  gr_mat_window_clear(H_top, ctx);
  gr_mat_window_clear(G_bottom, ctx);
  gr_mat_window_clear(H_bottom, ctx);
  return status;
}
