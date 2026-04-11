#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"

#include "addition.h"
#include "compression.h"
#include "flint/gr_types.h"
#include "matrix_aux.h"
#include "multiplication.h"

int gr_mat_split_quadrants(gr_mat_t G_a, gr_mat_t H_a, gr_mat_t G_b, gr_mat_t H_b, gr_mat_t G_c, gr_mat_t H_c,
                           gr_mat_t G_d, gr_mat_t H_d, gr_mat_t G_A, gr_mat_t H_A, gr_ctx_t ctx) {

  int status = GR_SUCCESS;
  slong n = gr_mat_nrows(G_A, ctx);
  slong rank = gr_mat_ncols(G_A, ctx);
  slong n1 = (n + 1) / 2;
  slong n2 = n / 2;

  gr_mat_init(G_a, n1, rank, ctx);
  gr_mat_init(H_a, n1, rank, ctx);

  gr_mat_t G_bottom, H_bottom;
  gr_mat_init(G_bottom, n2, rank, ctx);
  gr_mat_init(H_bottom, n2, rank, ctx);

  // horizontal slicing
  for (slong r = 0; r < n1; r++) {
    for (slong c = 0; c < rank; c++) {
      status |= gr_set(gr_mat_entry_ptr(G_a, r, c, ctx), gr_mat_entry_srcptr(G_A, r, c, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_a, r, c, ctx), gr_mat_entry_srcptr(H_A, r, c, ctx), ctx);
    }
  }
  for (slong r = 0; r < n2; r++) {
    for (slong c = 0; c < rank; c++) {
      status |= gr_set(gr_mat_entry_ptr(G_bottom, r, c, ctx), gr_mat_entry_srcptr(G_A, n1 + r, c, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_bottom, r, c, ctx), gr_mat_entry_srcptr(H_A, n1 + r, c, ctx), ctx);
    }
  }

  // e vectors
  gr_mat_t e_last, Za_elast, ZaT_elast, e0_n2_Hb, e0_n2_Gc;
  gr_mat_init(e_last, n1, 1, ctx);
  status |= gr_mat_zero(e_last, ctx);
  status |= gr_one(gr_mat_entry_ptr(e_last, n1 - 1, 0, ctx), ctx);

  gr_mat_init(Za_elast, n1, 1, ctx);
  gr_mat_init(ZaT_elast, n1, 1, ctx);

  {
    gr_mat_t tmp;
    gr_mat_init(tmp, n1, 1, ctx);
    // \Phi_+(b) correction: Z * a * e_{last}
    status |= gr_mat_mul_vector(tmp, G_a, H_a, e_last, ctx);
    status |= gr_mat_apply_Z(Za_elast, tmp, ctx);

    // \Phi_+(c) correction: Z * a^T * e_{last}
    status |= gr_mat_mul_vector(tmp, H_a, G_a, e_last, ctx);
    status |= gr_mat_apply_Z(ZaT_elast, tmp, ctx);
    gr_mat_clear(tmp, ctx);
  }

  gr_mat_init(e0_n2_Hb, n2, 1, ctx);
  status |= gr_mat_zero(e0_n2_Hb, ctx);
  status |= gr_one(gr_mat_entry_ptr(e0_n2_Hb, 0, 0, ctx), ctx);

  gr_mat_init(e0_n2_Gc, n2, 1, ctx);
  status |= gr_mat_zero(e0_n2_Gc, ctx);
  status |= gr_one(gr_mat_entry_ptr(e0_n2_Gc, 0, 0, ctx), ctx);

  // block b: G_top * H_bottom^T + (Z * a * e_{last}) * e_0^T
  slong rank_b = rank + 1;
  gr_mat_init(G_b, n1, rank_b, ctx);
  gr_mat_init(H_b, n2, rank_b, ctx);
  status |= gr_mat_concat_horizontal(G_b, G_a, Za_elast, ctx);
  status |= gr_mat_concat_horizontal(H_b, H_bottom, e0_n2_Hb, ctx);

  // block c: G_bottom * H_top^T + e_0 * (Z * a^T * e_{last})^T
  slong rank_c = rank + 1;
  gr_mat_init(G_c, n2, rank_c, ctx);
  gr_mat_init(H_c, n1, rank_c, ctx);
  status |= gr_mat_concat_horizontal(G_c, G_bottom, e0_n2_Gc, ctx);
  status |= gr_mat_concat_horizontal(H_c, H_a, ZaT_elast, ctx);

  // block: d
  gr_mat_t G_top_trunc, H_top_trunc;
  gr_mat_init(G_top_trunc, n2, rank, ctx);
  gr_mat_init(H_top_trunc, n2, rank, ctx);
  for (slong r = 0; r < n2; r++) {
    for (slong c = 0; c < rank; c++) {
      status |= gr_set(gr_mat_entry_ptr(G_top_trunc, r, c, ctx), gr_mat_entry_srcptr(G_a, r, c, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_top_trunc, r, c, ctx), gr_mat_entry_srcptr(H_a, r, c, ctx), ctx);
    }
  }

  status |= gr_mat_addition_generateur(G_bottom, H_bottom, G_top_trunc, H_top_trunc, G_d, H_d, ctx);

  gr_mat_clear(G_bottom, ctx);
  gr_mat_clear(H_bottom, ctx);
  gr_mat_clear(e_last, ctx);
  gr_mat_clear(Za_elast, ctx);
  gr_mat_clear(ZaT_elast, ctx);
  gr_mat_clear(e0_n2_Hb, ctx);
  gr_mat_clear(e0_n2_Gc, ctx);
  gr_mat_clear(G_top_trunc, ctx);
  gr_mat_clear(H_top_trunc, ctx);
  return status;
}



int gr_mat_pack_quadrants(gr_mat_t G_D, gr_mat_t H_D, gr_mat_t G_x, gr_mat_t H_x, gr_mat_t G_y, gr_mat_t H_y,
                          gr_mat_t G_z, gr_mat_t H_z, gr_mat_t G_t, gr_mat_t H_t, gr_ctx_t ctx) {

  int status = GR_SUCCESS;
  slong n1 = gr_mat_nrows(G_x, ctx);
  slong n2 = gr_mat_nrows(G_t, ctx);
  slong n = n1 + n2;

  slong rx = gr_mat_ncols(G_x, ctx);
  slong ry = gr_mat_ncols(G_y, ctx);
  slong rz = gr_mat_ncols(G_z, ctx);
  slong rt = gr_mat_ncols(G_t, ctx);
  slong total_rank = rx + ry + rz + rt;

  gr_mat_init(G_D, n, total_rank, ctx);
  gr_mat_init(H_D, n, total_rank, ctx);
  status |= gr_mat_zero(G_D, ctx);
  status |= gr_mat_zero(H_D, ctx);

  // not that elegantly pack the quadrants
  const gr_mat_struct *G_blk[4] = {G_x, G_y, G_z, G_t};
  const gr_mat_struct *H_blk[4] = {H_x, H_y, H_z, H_t};
  const slong G_rs[4] = {0, 0, n1, n1};
  const slong G_rows[4] = {n1, n1, n2, n2};
  const slong H_rs[4] = {0, n1, 0, n1};
  const slong H_rows[4] = {n1, n2, n1, n2};
  const slong col_off[4] = {0, rx, rx + ry, rx + ry + rz};

  for (slong b = 0; b < 4; b++) {
    slong rk = gr_mat_ncols(G_blk[b], ctx);
    for (slong k = 0; k < rk; k++) {
      slong col = col_off[b] + k;
      for (slong i = 0; i < G_rows[b]; i++) {
        status |= gr_set(gr_mat_entry_ptr(G_D, G_rs[b] + i, col, ctx), gr_mat_entry_srcptr(G_blk[b], i, k, ctx), ctx);
      }
      for (slong j = 0; j < H_rows[b]; j++) {
        status |= gr_set(gr_mat_entry_ptr(H_D, H_rs[b] + j, col, ctx), gr_mat_entry_srcptr(H_blk[b], j, k, ctx), ctx);
      }
    }
  }

  // reverse bleed
  gr_mat_t v_x, r_x, v_z, r_y, e_last_x, tmp_n1, tmp_n2;
  gr_mat_init(v_x, n1, 1, ctx);
  gr_mat_init(r_x, n1, 1, ctx);
  gr_mat_init(v_z, n2, 1, ctx);
  gr_mat_init(r_y, n2, 1, ctx);
  gr_mat_init(tmp_n1, n1, 1, ctx);
  gr_mat_init(tmp_n2, n2, 1, ctx);

  gr_mat_init(e_last_x, n1, 1, ctx);
  status |= gr_mat_zero(e_last_x, ctx);
  status |= gr_one(gr_mat_entry_ptr(e_last_x, n1 - 1, 0, ctx), ctx);

  gr_ptr s_val = gr_heap_init(ctx);

  // boundary bleed by Z
  status |= gr_mat_mul_vector(tmp_n1, G_x, H_x, e_last_x, ctx);
  status |= gr_set(s_val, gr_mat_entry_ptr(tmp_n1, n1 - 1, 0, ctx), ctx); // corner
  status |= gr_mat_apply_Z(v_x, tmp_n1, ctx);

  status |= gr_mat_mul_vector(tmp_n1, H_x, G_x, e_last_x, ctx);
  status |= gr_mat_apply_Z(r_x, tmp_n1, ctx);

  status |= gr_mat_mul_vector(tmp_n2, G_z, H_z, e_last_x, ctx);
  status |= gr_mat_apply_Z(v_z, tmp_n2, ctx);

  status |= gr_mat_mul_vector(tmp_n2, H_y, G_y, e_last_x, ctx);
  status |= gr_mat_apply_Z(r_y, tmp_n2, ctx);

  // cancel
  gr_mat_t G_corr, H_corr;
  gr_mat_init(G_corr, n, 2, ctx);
  gr_mat_init(H_corr, n, 2, ctx);
  status |= gr_mat_zero(G_corr, ctx);
  status |= gr_mat_zero(H_corr, ctx);

  gr_ptr temp_calc = gr_heap_init(ctx);

  for (slong i = 0; i < n1; i++) {
    // y bleed: -v_x * e0^T
    status |= gr_neg(temp_calc, gr_mat_entry_ptr(v_x, i, 0, ctx), ctx);
    status |= gr_set(gr_mat_entry_ptr(G_corr, i, 0, ctx), temp_calc, ctx);

    // z bleed: -e0 * r_x^T
    status |= gr_neg(temp_calc, gr_mat_entry_ptr(r_x, i, 0, ctx), ctx);
    status |= gr_set(gr_mat_entry_ptr(H_corr, i, 1, ctx), temp_calc, ctx);
  }

  for (slong i = 0; i < n2; i++) {
    // t bleed (from z and corner): (-v_z - s_val) * e0^T
    status |= gr_neg(temp_calc, gr_mat_entry_ptr(v_z, i, 0, ctx), ctx);
    if (i == 0) status |= gr_sub(temp_calc, temp_calc, s_val, ctx);
    status |= gr_set(gr_mat_entry_ptr(G_corr, n1 + i, 0, ctx), temp_calc, ctx);

    if (i == 0) {
      status |= gr_one(gr_mat_entry_ptr(H_corr, n1 + i, 0, ctx), ctx);
      status |= gr_one(gr_mat_entry_ptr(G_corr, n1 + i, 1, ctx), ctx);
    }

    // t bleed (from y): -e0 * r_y^T
    status |= gr_neg(temp_calc, gr_mat_entry_ptr(r_y, i, 0, ctx), ctx);
    status |= gr_set(gr_mat_entry_ptr(H_corr, n1 + i, 1, ctx), temp_calc, ctx);
  }

  gr_mat_t G_final, H_final;
  status |= gr_mat_addition_generateur(G_D, H_D, G_corr, H_corr, G_final, H_final, ctx);
  gr_mat_swap(G_D, G_final, ctx);
  gr_mat_swap(H_D, H_final, ctx);
  status |= gr_mat_generator_compress(G_D, H_D, ctx);

  gr_mat_clear(G_final, ctx);
  gr_mat_clear(H_final, ctx);
  gr_mat_clear(G_corr, ctx);
  gr_mat_clear(H_corr, ctx);
  gr_mat_clear(v_x, ctx);
  gr_mat_clear(r_x, ctx);
  gr_mat_clear(v_z, ctx);
  gr_mat_clear(r_y, ctx);
  gr_mat_clear(tmp_n1, ctx);
  gr_mat_clear(tmp_n2, ctx);
  gr_mat_clear(e_last_x, ctx);
  gr_heap_clear(s_val, ctx);
  gr_heap_clear(temp_calc, ctx);

  return status;
}