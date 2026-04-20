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

  slong total_rank = rx + ry + rz + rt + 2;

  gr_mat_init(G_D, n, total_rank, ctx);
  gr_mat_init(H_D, n, total_rank, ctx);
  status |= gr_mat_zero(G_D, ctx);
  status |= gr_mat_zero(H_D, ctx);

  slong col_x = 0;
  slong col_y = rx;
  slong col_z = rx + ry;
  slong col_t = rx + ry + rz;

  // X
  for (slong j = 0; j < rx; j++) {
    for (slong i = 0; i < n1; i++) {
      status |= gr_set(gr_mat_entry_ptr(G_D, i, col_x + j, ctx), gr_mat_entry_srcptr(G_x, i, j, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_D, i, col_x + j, ctx), gr_mat_entry_srcptr(H_x, i, j, ctx), ctx);
    }
  }
  // Y
  for (slong j = 0; j < ry; j++) {
    for (slong i = 0; i < n1; i++)
      status |= gr_set(gr_mat_entry_ptr(G_D, i, col_y + j, ctx), gr_mat_entry_srcptr(G_y, i, j, ctx), ctx);
    for (slong i = 0; i < n2; i++)
      status |= gr_set(gr_mat_entry_ptr(H_D, n1 + i, col_y + j, ctx), gr_mat_entry_srcptr(H_y, i, j, ctx), ctx);
  }
  // Z
  for (slong j = 0; j < rz; j++) {
    for (slong i = 0; i < n2; i++)
      status |= gr_set(gr_mat_entry_ptr(G_D, n1 + i, col_z + j, ctx), gr_mat_entry_srcptr(G_z, i, j, ctx), ctx);
    for (slong i = 0; i < n1; i++)
      status |= gr_set(gr_mat_entry_ptr(H_D, i, col_z + j, ctx), gr_mat_entry_srcptr(H_z, i, j, ctx), ctx);
  }
  // T
  for (slong j = 0; j < rt; j++) {
    for (slong i = 0; i < n2; i++) {
      status |= gr_set(gr_mat_entry_ptr(G_D, n1 + i, col_t + j, ctx), gr_mat_entry_srcptr(G_t, i, j, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_D, n1 + i, col_t + j, ctx), gr_mat_entry_srcptr(H_t, i, j, ctx), ctx);
    }
  }
  
  status |= gr_mat_generator_compress(G_D, H_D, ctx);
  
  return status;
}