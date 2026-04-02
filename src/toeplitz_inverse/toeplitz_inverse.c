#include "addition.h"
#include "compression.h"
#include "matrix_aux.h"
#include "multiplication.h"

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"

int gr_toeplitz_inverse(gr_mat_t G_D, gr_mat_t H_D, gr_mat_t G_A, gr_mat_t H_A, gr_ctx_t ctx) {

  /*
   * Algorithm 10.1 Algorithm of type Strassen to inverse a matrix
   * quasi-Toeplitz, from [1].
   *
   * All generators are the phi+ displacement operator:
   *   phi+(A) = A - Z*A*Z^T = G * H^T
   * We trust that at the end of the recursion phi+(A) = phi-(A) where A is a 1x1 matrix
   *
   * Block decomposition of A into [a b; c d] with generators (G_top, H_top)
   * for a, and the b/c/d generators derived from the report (section 3.3.3):
   *
   *   a : Ga = Gtop,                      Ha = Htop
   *   b : Gb = [Gtop | Z*a*e_{n1-1}],     Hb = [Hbottom | e0_{n1}]
   *   c : Gc = [Gbottom | e0_{n2}],       Hc = [Htop | Z*a^T*e_{n1-1}]
   *   d : Gd = [Gbottom | Gtop],          Hd = [Hbottom | Htop]
   *
   * where:
   *   e_{n1-1} is the last standard basis vector of size n1
   *   e0_{n1}  is the first standard basis vector of size n1  (for Hb)
   *   e0_{n2}  is the first standard basis vector of size n2  (for Gc)
   *   Z*a*e_{n1-1}   = apply_Z(mul_vector(G_top, H_top, e_{n1-1}))
   *   Z*a^T*e_{n1-1} = apply_Z(mul_vector(H_top, G_top, e_{n1-1}))
   *                    (a^T has generators (H_top, G_top))
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

  // 2.0: prepare to calculate the generators a b c d for n/2
  slong n1 = (n + 1) / 2;
  slong n2 = n / 2;

  gr_mat_t G_top, H_top, G_bottom, H_bottom;
  // gr_mat_window_init(G_top, G_A, 0, 0, n1, rank, ctx);
  // gr_mat_window_init(H_top, H_A, 0, 0, n1, rank, ctx);
  // gr_mat_window_init(G_bottom, G_A, n1, 0, n, rank, ctx);
  // gr_mat_window_init(H_bottom, H_A, n1, 0, n, rank, ctx);
  gr_mat_init(G_top, n1, rank, ctx);
  gr_mat_init(H_top, n1, rank, ctx);
  gr_mat_init(G_bottom, n2, rank, ctx);
  gr_mat_init(H_bottom, n2, rank, ctx);
  for (slong r = 0; r < n1; r++) {
    for (slong c = 0; c < rank; c++) {
      status |= gr_set(gr_mat_entry_ptr(G_top, r, c, ctx), gr_mat_entry_srcptr(G_A, r, c, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_top, r, c, ctx), gr_mat_entry_srcptr(H_A, r, c, ctx), ctx);
    }
  }
  for (slong r = 0; r < n2; r++) {
    for (slong c = 0; c < rank; c++) {
      status |= gr_set(gr_mat_entry_ptr(G_bottom, r, c, ctx), gr_mat_entry_srcptr(G_A, n1 + r, c, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_bottom, r, c, ctx), gr_mat_entry_srcptr(H_A, n1 + r, c, ctx), ctx);
    }
  }

  // 2.1: build correction vectors for b and c

  // e_{n1-1}: last basis vector of size n1, shared by both b and c
  gr_mat_t e_last, Za_elast, ZaT_elast, e0_n1, e0_n2;
  gr_mat_init(e_last, n1, 1, ctx);
  status |= gr_mat_zero(e_last, ctx);
  status |= gr_one(gr_mat_entry_ptr(e_last, n1 - 1, 0, ctx), ctx);

  // Z*a*e_{n1-1} — extra column for G_b
  gr_mat_init(Za_elast, n1, 1, ctx);
  {
    gr_mat_t tmp;
    gr_mat_init(tmp, n1, 1, ctx);
    status |= gr_mat_mul_vector(tmp, G_top, H_top, e_last, ctx);
    status |= gr_mat_apply_Z(Za_elast, tmp, ctx);
    gr_mat_clear(tmp, ctx);
  }

  // Z*a^T*e_{n1-1} — extra column for H_c (a^T = H_top, G_top)
  gr_mat_init(ZaT_elast, n1, 1, ctx);
  {
    gr_mat_t tmp;
    gr_mat_init(tmp, n1, 1, ctx);
    status |= gr_mat_mul_vector(tmp, H_top, G_top, e_last, ctx);
    status |= gr_mat_apply_Z(ZaT_elast, tmp, ctx);
    gr_mat_clear(tmp, ctx);
  }

  // e0 of size n1 — appended to H_b
  gr_mat_init(e0_n1, n1, 1, ctx);
  status |= gr_mat_zero(e0_n1, ctx);
  status |= gr_one(gr_mat_entry_ptr(e0_n1, 0, 0, ctx), ctx);

  // e0 of size n2 — appended to G_c
  gr_mat_init(e0_n2, n2, 1, ctx);
  status |= gr_mat_zero(e0_n2, ctx);
  status |= gr_one(gr_mat_entry_ptr(e0_n2, 0, 0, ctx), ctx);

  if (status != GR_SUCCESS) {
    gr_mat_clear(e_last, ctx);
    gr_mat_clear(Za_elast, ctx);
    gr_mat_clear(ZaT_elast, ctx);
    gr_mat_clear(e0_n1, ctx);
    gr_mat_clear(e0_n2, ctx);
    goto free_windows;
  }

  // 2.2: b : Gb = [Gtop | Za_elast],  Hb = [Hbottom | e0_n1]
  gr_mat_t G_b, H_b;
  slong rank_b = gr_mat_ncols(G_top, ctx) + gr_mat_ncols(Za_elast, ctx);
  gr_mat_init(G_b, n1, rank_b, ctx);
  gr_mat_init(H_b, n1, rank_b, ctx); // note: both have n1 rows
  status |= gr_mat_concat_horizontal(G_b, G_top, Za_elast, ctx);
  status |= gr_mat_concat_horizontal(H_b, H_bottom, e0_n1, ctx);
  gr_mat_clear(Za_elast, ctx);
  gr_mat_clear(e0_n1, ctx);

  // 2.3: c : Gc = [Gbottom | e0_n2],  Hc = [Htop | ZaT_elast]
  gr_mat_t G_c, H_c;
  slong rank_c = gr_mat_ncols(G_bottom, ctx) + gr_mat_ncols(e0_n2, ctx);
  gr_mat_init(G_c, n2, rank_c, ctx);
  gr_mat_init(H_c, n2, rank_c, ctx);
  status |= gr_mat_concat_horizontal(G_c, G_bottom, e0_n2, ctx);
  status |= gr_mat_concat_horizontal(H_c, H_top, ZaT_elast, ctx);
  gr_mat_clear(ZaT_elast, ctx);
  gr_mat_clear(e0_n2, ctx);
  gr_mat_clear(e_last, ctx);

  // 2.4: d : Gd = [Gbottom | Gtop],  Hd = [Hbottom | Htop]
  gr_mat_t G_d, H_d;
  status |= gr_mat_addition_generateur(G_bottom, H_bottom, G_top, H_top, G_d, H_d, ctx);
  if (status != GR_SUCCESS) {
    gr_mat_clear(G_d, ctx);
    gr_mat_clear(H_d, ctx);
    gr_mat_clear(G_c, ctx);
    gr_mat_clear(H_c, ctx);
    gr_mat_clear(G_b, ctx);
    gr_mat_clear(H_b, ctx);
    goto free_windows;
  }

  // 3: calculate rec e := a^{-1}
  gr_mat_t G_e, H_e;
  gr_mat_init(G_e, n1, rank, ctx);
  gr_mat_init(H_e, n1, rank, ctx);
  status |= gr_toeplitz_inverse(G_e, H_e, G_top, H_top, ctx);
  if (status != GR_SUCCESS) goto free_e_top_bottom;

  // 4: calculate generators for S := d - c*e*b

  // 4.1 ce = c * e
  gr_mat_t G_ce, H_ce;
  status |= gr_mat_mul_generator(G_ce, H_ce, G_c, H_c, G_e, H_e, ctx);
  gr_mat_clear(G_c, ctx);
  gr_mat_clear(H_c, ctx);
  if (status != GR_SUCCESS) {
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_e_top_bottom;
  }
  status |= gr_mat_generator_compress(G_ce, H_ce, ctx);
  if (status != GR_SUCCESS) {
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_e_top_bottom;
  }

  // 4.2 eb = e * b
  gr_mat_t G_eb, H_eb;
  status |= gr_mat_mul_generator(G_eb, H_eb, G_e, H_e, G_b, H_b, ctx);
  gr_mat_clear(G_b, ctx);
  gr_mat_clear(H_b, ctx);
  if (status != GR_SUCCESS) {
    gr_mat_clear(G_eb, ctx);
    gr_mat_clear(H_eb, ctx);
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_e_top_bottom;
  }
  status |= gr_mat_generator_compress(G_eb, H_eb, ctx);
  if (status != GR_SUCCESS) {
    gr_mat_clear(G_eb, ctx);
    gr_mat_clear(H_eb, ctx);
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_e_top_bottom;
  }

  // 4.3 ceb = ce * eb, then negate G in-place
  gr_mat_t G_ceb, H_ceb;
  status |= gr_mat_mul_generator(G_ceb, H_ceb, G_ce, H_ce, G_eb, H_eb, ctx);
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
  status |= gr_mat_generator_compress(G_ceb, H_ceb, ctx);
  if (status != GR_SUCCESS) {
    gr_mat_clear(G_ceb, ctx);
    gr_mat_clear(H_ceb, ctx);
    gr_mat_clear(G_eb, ctx);
    gr_mat_clear(H_eb, ctx);
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_e_top_bottom;
  }

  // 4.4 S = d + (-ceb)
  gr_mat_t G_S, H_S;
  status |= gr_mat_addition_generateur(G_d, H_d, G_ceb, H_ceb, G_S, H_S, ctx);
  gr_mat_clear(G_ceb, ctx);
  gr_mat_clear(H_ceb, ctx);
  gr_mat_clear(G_d, ctx);
  gr_mat_clear(H_d, ctx);
  if (status != GR_SUCCESS) {
    gr_mat_clear(G_S, ctx);
    gr_mat_clear(H_S, ctx);
    gr_mat_clear(G_eb, ctx);
    gr_mat_clear(H_eb, ctx);
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_e_top_bottom;
  }
  status |= gr_mat_generator_compress(G_S, H_S, ctx);
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

  // 6: return the inv generators for A^{-1} = [x y; z t]
  //    x = e + ebt*ce,  y = -eb*t,  z = -t*ce,  t = t

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
  status |= gr_mat_generator_compress(G_ebt, H_ebt, ctx);
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
  status |= gr_mat_generator_compress(G_tce, H_tce, ctx);
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
  status |= gr_mat_generator_compress(G_ebtce, H_ebtce, ctx);
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
  status |= gr_mat_generator_compress(G_x, H_x, ctx);
  if (status != GR_SUCCESS) { goto free_x; }

  // y = -ebt (negate G in-place)
  status |= gr_mat_neg(G_ebt, G_ebt, ctx);
  if (status != GR_SUCCESS) { goto free_x; }

  // z = -tce (negate G in-place)
  status |= gr_mat_neg(G_tce, G_tce, ctx);
  if (status != GR_SUCCESS) { goto free_x; }

  // pack [x y; z t] into G_D / H_D
  {
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
free_windows:
  gr_mat_clear(G_top, ctx);
  gr_mat_clear(H_top, ctx);
  gr_mat_clear(G_bottom, ctx);
  gr_mat_clear(H_bottom, ctx);
  return status;
}