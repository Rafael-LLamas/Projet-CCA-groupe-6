#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"

#include "compression.h"
#include "flint/gr_types.h"

int gr_mat_split_quadrants(gr_mat_t G_a, gr_mat_t H_a, gr_mat_t G_b, gr_mat_t H_b, gr_mat_t G_c, gr_mat_t H_c,
                           gr_mat_t G_d, gr_mat_t H_d, gr_mat_t G_A, gr_mat_t H_A, gr_ctx_t ctx) {

  int status = GR_SUCCESS;
  slong n = gr_mat_nrows(G_A, ctx);
  slong rank = gr_mat_ncols(G_A, ctx);
  slong n1 = (n + 1) / 2;
  slong n2 = n / 2;

  gr_mat_init(G_a, n1, rank, ctx); // Top-Left
  gr_mat_init(H_a, n1, rank, ctx);
  gr_mat_init(G_b, n1, rank, ctx); // Top-Right
  gr_mat_init(H_b, n2, rank, ctx);
  gr_mat_init(G_c, n2, rank, ctx); // Bottom-Left
  gr_mat_init(H_c, n1, rank, ctx);
  gr_mat_init(G_d, n2, rank, ctx); // Bottom-Right
  gr_mat_init(H_d, n2, rank, ctx);

  for (slong r = 0; r < n1; r++) {
    for (slong c = 0; c < rank; c++) {
      status |= gr_set(gr_mat_entry_ptr(G_a, r, c, ctx), gr_mat_entry_srcptr(G_A, r, c, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_a, r, c, ctx), gr_mat_entry_srcptr(H_A, r, c, ctx), ctx);

      status |= gr_set(gr_mat_entry_ptr(G_b, r, c, ctx), gr_mat_entry_srcptr(G_A, r, c, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_c, r, c, ctx), gr_mat_entry_srcptr(H_A, r, c, ctx), ctx);
    }
  }

  for (slong r = 0; r < n2; r++) {
    for (slong c = 0; c < rank; c++) {
      status |= gr_set(gr_mat_entry_ptr(H_b, r, c, ctx), gr_mat_entry_srcptr(H_A, n1 + r, c, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(G_c, r, c, ctx), gr_mat_entry_srcptr(G_A, n1 + r, c, ctx), ctx);

      status |= gr_set(gr_mat_entry_ptr(G_d, r, c, ctx), gr_mat_entry_srcptr(G_A, n1 + r, c, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_d, r, c, ctx), gr_mat_entry_srcptr(H_A, n1 + r, c, ctx), ctx);
    }
  }

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

  slong col_x = 0;
  slong col_y = rx;
  slong col_z = rx + ry;
  slong col_t = rx + ry + rz;

  // x
  for (slong j = 0; j < rx; j++) {
    for (slong i = 0; i < n1; i++) {
      status |= gr_set(gr_mat_entry_ptr(G_D, i, col_x + j, ctx), gr_mat_entry_srcptr(G_x, i, j, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_D, i, col_x + j, ctx), gr_mat_entry_srcptr(H_x, i, j, ctx), ctx);
    }
  }
  // y
  for (slong j = 0; j < ry; j++) {
    for (slong i = 0; i < n1; i++)
      status |= gr_set(gr_mat_entry_ptr(G_D, i, col_y + j, ctx), gr_mat_entry_srcptr(G_y, i, j, ctx), ctx);
    for (slong i = 0; i < n2; i++)
      status |= gr_set(gr_mat_entry_ptr(H_D, n1 + i, col_y + j, ctx), gr_mat_entry_srcptr(H_y, i, j, ctx), ctx);
  }
  // z
  for (slong j = 0; j < rz; j++) {
    for (slong i = 0; i < n2; i++)
      status |= gr_set(gr_mat_entry_ptr(G_D, n1 + i, col_z + j, ctx), gr_mat_entry_srcptr(G_z, i, j, ctx), ctx);
    for (slong i = 0; i < n1; i++)
      status |= gr_set(gr_mat_entry_ptr(H_D, i, col_z + j, ctx), gr_mat_entry_srcptr(H_z, i, j, ctx), ctx);
  }
  // t
  for (slong j = 0; j < rt; j++) {
    for (slong i = 0; i < n2; i++) {
      status |= gr_set(gr_mat_entry_ptr(G_D, n1 + i, col_t + j, ctx), gr_mat_entry_srcptr(G_t, i, j, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_D, n1 + i, col_t + j, ctx), gr_mat_entry_srcptr(H_t, i, j, ctx), ctx);
    }
  }

  status |= gr_mat_generator_compress(G_D, H_D, ctx);
  
  return status;
}