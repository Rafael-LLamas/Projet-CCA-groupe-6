#include "compression.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"

// HELPERS ---------------------------------

// shift vector iterative (opti)
int shift_vec(gr_mat_t out, gr_mat_t v, slong n, gr_ctx_t ctx) {
  int status = GR_SUCCESS;
  status |= gr_mat_zero(out, ctx);
  for (slong i = 1; i < n; i++)
      
    status |= gr_set(gr_mat_entry_ptr(out, i, 0, ctx), gr_mat_entry_srcptr(v, i - 1, 0, ctx), ctx);
  return status;
}

// reconstructs a single cell (i, j) from generators G and H
int get_disp_cell(gr_ptr res, gr_mat_t G, gr_mat_t H, slong i, slong j, gr_ctx_t ctx) {
  int status = GR_SUCCESS;
  slong rank = gr_mat_ncols(G, ctx);
  slong limit = (i < j) ? i : j;

  gr_ptr tmp = gr_heap_init(ctx);
  status |= gr_zero(res, ctx);

  for (slong k = 0; k < rank; k++) {
    for (slong x = 0; x <= limit; x++) {
      status |= gr_mul(tmp, gr_mat_entry_srcptr(G, i - x, k, ctx), gr_mat_entry_srcptr(H, j - x, k, ctx), ctx);
      status |= gr_add(res, res, tmp, ctx);
    }
  }

  gr_heap_clear(tmp, ctx);
  return status;
}

// calculate v_a
int calculate_v_a(gr_mat_t v_a, gr_mat_t G, gr_mat_t H, slong n1, slong rank, gr_ctx_t ctx) {
  int status = GR_SUCCESS;
  for (slong i = 0; i < n1; i++) status |= get_disp_cell(gr_mat_entry_ptr(v_a, i, 0, ctx), G, H, i, n1 - 1, ctx);
  return status;
}

// calculate r_a
int calculate_r_a(gr_mat_t r_a, gr_mat_t G, gr_mat_t H, slong n1, slong rank, gr_ctx_t ctx) {
  int status = GR_SUCCESS;
  for (slong j = 0; j < n1; j++) status |= get_disp_cell(gr_mat_entry_ptr(r_a, j, 0, ctx), G, H, n1 - 1, j, ctx);
  return status;
}

// calculate scalar s_a
int calculate_s_a(gr_ptr s_a, gr_mat_t G, gr_mat_t H, slong n1, slong rank, gr_ctx_t ctx) {
  return get_disp_cell(s_a, G, H, n1 - 1, n1 - 1, ctx);
}

// calculate v_c (coords: n1 + l)
int calculate_v_c(gr_mat_t v_c, gr_mat_t G, gr_mat_t H, slong n1, slong n2, slong rank, gr_ctx_t ctx) {
  int status = GR_SUCCESS;
  for (slong l = 0; l < n2; l++) status |= get_disp_cell(gr_mat_entry_ptr(v_c, l, 0, ctx), G, H, n1 + l, n1 - 1, ctx);
  return status;
}

// calculate r_b (coords: n1 + m)
int calculate_r_b(gr_mat_t r_b, gr_mat_t G, gr_mat_t H, slong n1, slong n2, slong rank, gr_ctx_t ctx) {
  int status = GR_SUCCESS;
  for (slong m = 0; m < n2; m++) status |= get_disp_cell(gr_mat_entry_ptr(r_b, m, 0, ctx), G, H, n1 - 1, n1 + m, ctx);
  return status;
}

// MAIN SPLIT ---------------------------------

int gr_mat_split_quadrants(gr_mat_t G_a, gr_mat_t H_a, gr_mat_t G_b, gr_mat_t H_b, gr_mat_t G_c, gr_mat_t H_c,
                           gr_mat_t G_d, gr_mat_t H_d, gr_mat_t G_A, gr_mat_t H_A, gr_ctx_t ctx) {
  int status = GR_SUCCESS;

  slong n = gr_mat_nrows(G_A, ctx);
  slong rank = gr_mat_ncols(G_A, ctx);
  slong n1 = (n + 1) / 2;
  slong n2 = n / 2;

  // --- UTILITIES ---
  // vectors & scalar
  gr_mat_t v_a, r_a, v_c, r_b;
  gr_mat_init(v_a, n1, 1, ctx);
  gr_mat_init(r_a, n1, 1, ctx);
  gr_mat_init(v_c, n2, 1, ctx);
  gr_mat_init(r_b, n2, 1, ctx);
  status |= calculate_v_a(v_a, G_A, H_A, n1, rank, ctx);
  status |= calculate_r_a(r_a, G_A, H_A, n1, rank, ctx);
  status |= calculate_v_c(v_c, G_A, H_A, n1, n2, rank, ctx);
  status |= calculate_r_b(r_b, G_A, H_A, n1, n2, rank, ctx);

  gr_ptr s_a = gr_heap_init(ctx);
  status |= calculate_s_a(s_a, G_A, H_A, n1, rank, ctx);

  // shifted vectors
  gr_mat_t Zva, Zra, Zvc, Zrb;
  gr_mat_init(Zva, n1, 1, ctx);
  gr_mat_init(Zra, n1, 1, ctx);
  gr_mat_init(Zvc, n2, 1, ctx);
  gr_mat_init(Zrb, n2, 1, ctx);
  status |= shift_vec(Zva, v_a, n1, ctx);
  status |= shift_vec(Zra, r_a, n1, ctx);
  status |= shift_vec(Zvc, v_c, n2, ctx);
  status |= shift_vec(Zrb, r_b, n2, ctx);

  // e vectors
  gr_mat_t e_0_n1, e_0_n2;
  gr_mat_init(e_0_n1, n1, 1, ctx);
  gr_mat_init(e_0_n2, n2, 1, ctx);
  status |= gr_mat_zero(e_0_n1, ctx);
  status |= gr_mat_zero(e_0_n2, ctx);
  status |= gr_one(gr_mat_entry_ptr(e_0_n1, 0, 0, ctx), ctx);
  status |= gr_one(gr_mat_entry_ptr(e_0_n2, 0, 0, ctx), ctx);

  // G_bottom and H_bottom
  gr_mat_t G_bot, H_bot;
  gr_mat_init(G_bot, n2, rank, ctx);
  gr_mat_init(H_bot, n2, rank, ctx);
  for (slong r = 0; r < n2; r++) {
    for (slong c = 0; c < rank; c++) {
      status |= gr_set(gr_mat_entry_ptr(G_bot, r, c, ctx), gr_mat_entry_srcptr(G_A, n1 + r, c, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_bot, r, c, ctx), gr_mat_entry_srcptr(H_A, n1 + r, c, ctx), ctx);
    }
  }

  // --- BLOCKS ---
  // Block a : G_a = G_TOP,  H_a = H_TOP
  gr_mat_init(G_a, n1, rank, ctx);
  gr_mat_init(H_a, n1, rank, ctx);
  for (slong r = 0; r < n1; r++) {
    for (slong c = 0; c < rank; c++) {
      status |= gr_set(gr_mat_entry_ptr(G_a, r, c, ctx), gr_mat_entry_srcptr(G_A, r, c, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_a, r, c, ctx), gr_mat_entry_srcptr(H_A, r, c, ctx), ctx);
    }
  }

  // Block b : G_b = [G_a | Z*v_a],  H_b = [H_bot | e_0_n2]
  gr_mat_init(G_b, n1, rank + 1, ctx);
  gr_mat_init(H_b, n2, rank + 1, ctx);
  status |= gr_mat_concat_horizontal(G_b, G_a, Zva, ctx);
  status |= gr_mat_concat_horizontal(H_b, H_bot, e_0_n2, ctx);

  // Block c : G_c = [G_bot | e_0_n2],  H_c = [H_a | Z*r_a]
  gr_mat_init(G_c, n2, rank + 1, ctx);
  gr_mat_init(H_c, n1, rank + 1, ctx);
  status |= gr_mat_concat_horizontal(G_c, G_bot, e_0_n2, ctx);
  status |= gr_mat_concat_horizontal(H_c, H_a, Zra, ctx);

  // Block d : G_d = [G_BOTTOM | s_a*e_0_n2 + Z*v_c | e_0_n2],  H_d = [H_BOTTOM | e_0_n2 | Z*r_b]
  gr_mat_t M; // M = s_a*e_0_n2 + Z*v_c
  gr_mat_init(M, n2, 1, ctx);
  for (slong i = 0; i < n2; i++)
    status |= gr_set(gr_mat_entry_ptr(M, i, 0, ctx), gr_mat_entry_srcptr(Zvc, i, 0, ctx), ctx);
  status |= gr_add(gr_mat_entry_ptr(M, 0, 0, ctx), gr_mat_entry_srcptr(M, 0, 0, ctx), s_a, ctx);

  gr_mat_t G_d_tmp, H_d_tmp;
  gr_mat_init(G_d_tmp, n2, rank + 1, ctx);
  gr_mat_init(H_d_tmp, n2, rank + 1, ctx);

  // G_d = [G_BOTTOM | M | e_0_n2]
  gr_mat_init(G_d, n2, rank + 2, ctx);
  status |= gr_mat_concat_horizontal(G_d_tmp, G_bot, M, ctx);
  status |= gr_mat_concat_horizontal(G_d, G_d_tmp, e_0_n2, ctx);

  // H_d = [H_BOTTOM | e_0_n2 | Z*r_b]
  gr_mat_init(H_d, n2, rank + 2, ctx);
  status |= gr_mat_concat_horizontal(H_d_tmp, H_bot, e_0_n2, ctx);
  status |= gr_mat_concat_horizontal(H_d, H_d_tmp, Zrb, ctx);

  gr_mat_clear(v_a, ctx);
  gr_mat_clear(r_a, ctx);
  gr_mat_clear(v_c, ctx);
  gr_mat_clear(r_b, ctx);
  gr_mat_clear(Zva, ctx);
  gr_mat_clear(Zra, ctx);
  gr_mat_clear(Zvc, ctx);
  gr_mat_clear(Zrb, ctx);
  gr_mat_clear(e_0_n1, ctx);
  gr_mat_clear(e_0_n2, ctx);
  gr_mat_clear(G_bot, ctx);
  gr_mat_clear(H_bot, ctx);
  gr_mat_clear(M, ctx);
  gr_mat_clear(G_d_tmp, ctx);
  gr_mat_clear(H_d_tmp, ctx);
  gr_clear(s_a, ctx);
  return status;
}

// MAIN PACK ---------------------------------

/*
 * G_D:
 *         +-------+-------+-------+-------+--------+---------+--------------------+---------+
 *         |  rx   |  ry   |  rz   |  rt   |   c1   |   c2    |         c3         |   c4    |
 * +-------+-------+-------+-------+-------+--------+---------+--------------------+---------+
 * |       |       |       |       |       |        |         |                    |         |
 * | n1    |  G_x  |  G_y  |   0   |   0   |-Z*v_a  |    0    |         0          |    0    |
 * |       |       |       |       |       |        |         |                    |         |
 * +-------+-------+-------+-------+-------+--------+---------+--------------------+---------+
 * |       |       |       |       |       |        |         |                    |         |
 * | n2    |   0   |   0   |  G_z  |  G_t  |   0    |-e_0_n2  | -(s_a*e_0 + Z*v_c) |-e_0_n2  |
 * |       |       |       |       |       |        |         |                    |         |
 * +-------+-------+-------+-------+-------+--------+---------+--------------------+---------+
 *
 *
 * H_D:
 *       +-------+-------+-------+-------+--------+---------+--------------------+---------+
 *       |  rx   |  ry   |  rz   |  rt   |   c1   |   c2    |         c3         |   c4    |
 * +-----+-------+-------+-------+-------+--------+---------+--------------------+---------+
 * |     |       |       |       |       |        |         |                    |         |
 * | n1  |  H_x  |   0   |  H_z  |   0   |   0    |  Z*r_a  |         0          |    0    |
 * |     |       |       |       |       |        |         |                    |         |
 * +-----+-------+-------+-------+-------+--------+---------+--------------------+---------+
 * |     |       |       |       |       |        |         |                    |         |
 * | n2  |   0   |  H_y  |   0   |  H_t  | e_0_n2 |    0    |       e_0_n2       |  Z*r_b  |
 * |     |       |       |       |       |        |         |                    |         |
 * +-----+-------+-------+-------+-------+--------+---------+--------------------+---------+
 */

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

  slong total_rank = rx + ry + rz + rt + 4; // +4 for correction columns c

  gr_mat_init(G_D, n, total_rank, ctx);
  gr_mat_init(H_D, n, total_rank, ctx);
  status |= gr_mat_zero(G_D, ctx);
  status |= gr_mat_zero(H_D, ctx);

  slong col_x = 0;
  slong col_y = rx;
  slong col_z = rx + ry;
  slong col_t = rx + ry + rz;
  slong col_c = rx + ry + rz + rt; // correction columns starts

  // Block x: G_TOP, H_TOP
  for (slong j = 0; j < rx; j++) {
    for (slong i = 0; i < n1; i++) {
      status |= gr_set(gr_mat_entry_ptr(G_D, i, col_x + j, ctx), gr_mat_entry_srcptr(G_x, i, j, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_D, i, col_x + j, ctx), gr_mat_entry_srcptr(H_x, i, j, ctx), ctx);
    }
  }
  // Block y: G_TOP, H_BOTTOM
  for (slong j = 0; j < ry; j++) {
    for (slong i = 0; i < n1; i++)
      status |= gr_set(gr_mat_entry_ptr(G_D, i, col_y + j, ctx), gr_mat_entry_srcptr(G_y, i, j, ctx), ctx);
    for (slong i = 0; i < n2; i++)
      status |= gr_set(gr_mat_entry_ptr(H_D, n1 + i, col_y + j, ctx), gr_mat_entry_srcptr(H_y, i, j, ctx), ctx);
  }
  // Block z: G_BOTTOM, H_TOP
  for (slong j = 0; j < rz; j++) {
    for (slong i = 0; i < n2; i++)
      status |= gr_set(gr_mat_entry_ptr(G_D, n1 + i, col_z + j, ctx), gr_mat_entry_srcptr(G_z, i, j, ctx), ctx);
    for (slong i = 0; i < n1; i++)
      status |= gr_set(gr_mat_entry_ptr(H_D, i, col_z + j, ctx), gr_mat_entry_srcptr(H_z, i, j, ctx), ctx);
  }
  // Block t: G_BOTTOM, H_BOTTOM
  for (slong j = 0; j < rt; j++) {
    for (slong i = 0; i < n2; i++) {
      status |= gr_set(gr_mat_entry_ptr(G_D, n1 + i, col_t + j, ctx), gr_mat_entry_srcptr(G_t, i, j, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_D, n1 + i, col_t + j, ctx), gr_mat_entry_srcptr(H_t, i, j, ctx), ctx);
    }
  }

  // bleed vectors
  gr_mat_t v_a, r_a, v_c, r_b;
  gr_mat_init(v_a, n1, 1, ctx);
  gr_mat_init(r_a, n1, 1, ctx);
  gr_mat_init(v_c, n2, 1, ctx);
  gr_mat_init(r_b, n2, 1, ctx);
  gr_ptr s_a = gr_heap_init(ctx);

  status |= calculate_v_a(v_a, G_x, H_x, n1, rx, ctx); // v_a from block x
  status |= calculate_r_a(r_a, G_x, H_x, n1, rx, ctx); // r_a from block x
  status |= calculate_s_a(s_a, G_x, H_x, n1, rx, ctx); // s_a from block x

  for (slong l = 0; l < n2; l++) { // v_c from block z
    status |= get_disp_cell(gr_mat_entry_ptr(v_c, l, 0, ctx), G_z, H_z, l, n1 - 1, ctx);
  }

  for (slong m = 0; m < n2; m++) { // r_b from block y
    status |= get_disp_cell(gr_mat_entry_ptr(r_b, m, 0, ctx), G_y, H_y, n1 - 1, m, ctx);
  }

  // filling to G_D and H_D
  slong c1 = col_c;
  slong c2 = col_c + 1;
  slong c3 = col_c + 2;
  slong c4 = col_c + 3;

  gr_ptr minus_1 = gr_heap_init(ctx);
  status |= gr_one(minus_1, ctx);
  status |= gr_neg(minus_1, minus_1, ctx);

  gr_mat_t Zva, Zra, Zvc, Zrb;
  gr_mat_init(Zva, n1, 1, ctx);
  gr_mat_init(Zra, n1, 1, ctx);
  gr_mat_init(Zvc, n2, 1, ctx);
  gr_mat_init(Zrb, n2, 1, ctx);
  status |= shift_vec(Zva, v_a, n1, ctx);
  status |= shift_vec(Zra, r_a, n1, ctx);
  status |= shift_vec(Zvc, v_c, n2, ctx);
  status |= shift_vec(Zrb, r_b, n2, ctx);
  
  
  gr_ptr tmp = gr_heap_init(ctx);

  // c1: G=-Z*v_a (top), H=e_0_n2 (bottom)
  for (slong i = 0; i < n1; i++) {
    status |= gr_neg(tmp, gr_mat_entry_srcptr(Zva, i, 0, ctx), ctx);
    status |= gr_set(gr_mat_entry_ptr(G_D, i, c1, ctx), tmp, ctx);
  }
  status |= gr_one(gr_mat_entry_ptr(H_D, n1, c1, ctx), ctx); // e_0_n2

  // c2: G=-e_0_n2 (bottom), H=Z*r_a (top)
  status |= gr_set(gr_mat_entry_ptr(G_D, n1, c2, ctx), minus_1, ctx); // -e_0_n2
  for (slong i = 0; i < n1; i++) {
    status |= gr_set(gr_mat_entry_ptr(H_D, i, c2, ctx), gr_mat_entry_srcptr(Zra, i, 0, ctx), ctx);
  }

  // c3: G=-(s_a*e_0 + Z*v_c) (bottom), H=e_0_n2 (bottom)
  for (slong i = 0; i < n2; i++) {
    if (i == 0) {
      status |= gr_add(tmp, s_a, gr_mat_entry_srcptr(Zvc, 0, 0, ctx), ctx);
    } else {
      status |= gr_set(tmp, gr_mat_entry_srcptr(Zvc, i, 0, ctx), ctx);
    }
    status |= gr_neg(tmp, tmp, ctx);
    status |= gr_set(gr_mat_entry_ptr(G_D, n1 + i, c3, ctx), tmp, ctx);
  }
  status |= gr_one(gr_mat_entry_ptr(H_D, n1, c3, ctx), ctx); // e_0_n2

  // c4: G=-e_0_n2 (bottom), H=Z*r_b (bottom)
  status |= gr_set(gr_mat_entry_ptr(G_D, n1, c4, ctx), minus_1, ctx); // -e_0_n2
  for (slong i = 0; i < n2; i++) {
    status |= gr_set(gr_mat_entry_ptr(H_D, n1 + i, c4, ctx), gr_mat_entry_srcptr(Zrb, i, 0, ctx), ctx);
  }

  gr_mat_clear(v_a, ctx);
  gr_mat_clear(r_a, ctx);
  gr_mat_clear(v_c, ctx);
  gr_mat_clear(r_b, ctx);
  gr_heap_clear(s_a, ctx);
  gr_heap_clear(minus_1, ctx);
  gr_heap_clear(tmp, ctx);

  status |= gr_mat_generator_compress(G_D, H_D, ctx);

  return status;
}