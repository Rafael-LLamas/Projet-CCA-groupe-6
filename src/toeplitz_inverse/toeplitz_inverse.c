#include "addition.h"
#include "compression.h"
#include "matrix_aux.h"
#include "multiplication.h"

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"

int gr_toeplitz_inverse(gr_mat_t G_D, gr_mat_t H_D, gr_mat_t G_A, gr_mat_t H_A, gr_ctx_t ctx) {

  /*
   * The core idea of this algorithm is the Strassen Block Inversion,
   * instead of inverting a large matrix, I can cut down the matrix to
   * 4 pieces: [[a ,b] [c, d]].
   * The catch is that we never touch the dense matrices themselves, every
   * operation is made on their generators.
   */

  int status = GR_SUCCESS;
  slong n = gr_mat_nrows(G_A, ctx);
  slong rank = gr_mat_ncols(G_A, ctx);

  // TODO: can maybe add a check of its determinant once implemented the toeplitz version (det must not be zero)
  printf("[inv] n=%ld rank=%ld\n", n, rank);

  /* BASE CASE: -------------------------------------------
   * If a matrix is 1x1, its shift operator Z is 0.
   * Therefore its displacement DISP_PLUS and DISP_MINS is the same.
   *
   * Once recieved, we simply invert the value inside rthe G_A and
   * H_A and place it to G_D and H_D.
   */
  if (n == 1) {
    if (rank == 0) {
      printf("[inv] BASE CASE FAILED: rank=0\n");
      return GR_UNABLE;
    }
    printf("[inv] base case G_A[0,0]=");
    gr_print(gr_mat_entry_srcptr(G_A, 0, 0, ctx), ctx);
    printf(" H_A[0,0]=");
    gr_print(gr_mat_entry_srcptr(H_A, 0, 0, ctx), ctx);
    printf("\n");
    status |= gr_mat_zero(G_D, ctx);
    status |= gr_mat_zero(H_D, ctx);
    status |= gr_inv(gr_mat_entry_ptr(G_D, 0, 0, ctx), gr_mat_entry_ptr(G_A, 0, 0, ctx), ctx);
    status |= gr_inv(gr_mat_entry_ptr(H_D, 0, 0, ctx), gr_mat_entry_ptr(H_A, 0, 0, ctx), ctx);
    if (status != GR_SUCCESS) printf("[inv] BASE CASE: gr_inv failed\n");
    return status;
  }

  /*
   * SPLITTING THE MATRIX IN HALF: -------------------------------------------
   * We divide our nxn matrix to 4 quadrants: a,b,c,d.
   * To accomplish this division we slice G_A and H_A horizontally
   *
   * The top n1 corresponds to the top blocks a and d,
   * The bottom n2 corresponds to the bottom blocks b and c.
   *
   * We calculate n1 and n2 then allocate memory for G_top,
   * H_top, G_bottom, H_bottom.
   * Then copy respective rows from the input generators.
   */

  // 2.0: prepare to calculate the generators a b c d for n/2
  slong n1 = (n + 1) / 2;
  slong n2 = n / 2;
  printf("[inv] n1=%ld n2=%ld\n", n1, n2);

  gr_mat_t G_top, H_top, G_bottom, H_bottom;
  gr_mat_init(G_top, n1, rank, ctx);
  gr_mat_init(H_top, n1, rank, ctx);
  gr_mat_init(G_bottom, n2, rank, ctx);
  gr_mat_init(H_bottom, n2, rank, ctx);
  for (slong r = 0; r < n1; r++)
    for (slong c = 0; c < rank; c++) {
      status |= gr_set(gr_mat_entry_ptr(G_top, r, c, ctx), gr_mat_entry_srcptr(G_A, r, c, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_top, r, c, ctx), gr_mat_entry_srcptr(H_A, r, c, ctx), ctx);
    }
  for (slong r = 0; r < n2; r++)
    for (slong c = 0; c < rank; c++) {
      status |= gr_set(gr_mat_entry_ptr(G_bottom, r, c, ctx), gr_mat_entry_srcptr(G_A, n1 + r, c, ctx), ctx);
      status |= gr_set(gr_mat_entry_ptr(H_bottom, r, c, ctx), gr_mat_entry_srcptr(H_A, n1 + r, c, ctx), ctx);
    }

  if (status != GR_SUCCESS) {
    printf("[inv] FAILED: splitting G_top/H_top/G_bottom/H_bottom\n");
    goto free_windows;
  }

  printf("[inv] G_top:\n");
  gr_mat_print(G_top, ctx);
  printf("[inv] H_top:\n");
  gr_mat_print(H_top, ctx);
  printf("[inv] G_bottom:\n");
  gr_mat_print(G_bottom, ctx);
  printf("[inv] H_bottom:\n");
  gr_mat_print(H_bottom, ctx);

  /* DISPLACEMENT BLEED: -------------------------------------------
   * As we are working with generators, once a reconstruction
   * is called, we do diagonal operations to retrive the matrix from
   * its displacement. (As generators are formed from the displacement
   * matrix, the values are not from the original matrix)
   *
   * For DISP_PLUS, we shift the matrix one down and right.
   * The bottom right corner of a bleeds into d and edges bleeds
   * into b and c.
   *
   * Detailed math can be found in the report.
   * We create a standart basis vectors e_last and e_0.
   * We multiply them by our generators to extract specific columns and
   * rows we need and apply Z shift to simulate the diagonal bleed.
   */

  // 2.1: build correction vectors for b and c
  gr_mat_t e_last, Za_elast, ZaT_elast, e0_n2_Hb, e0_n2_Gc;

  // create a vector with a '1' at the very bottom (size n1)
  gr_mat_init(e_last, n1, 1, ctx);
  status |= gr_mat_zero(e_last, ctx);
  status |= gr_one(gr_mat_entry_ptr(e_last, n1 - 1, 0, ctx), ctx);

  // calculate Z * a * e_{last}  (bleed for b)
  gr_mat_init(Za_elast, n1, 1, ctx);
  {
    gr_mat_t tmp;
    gr_mat_init(tmp, n1, 1, ctx);
    status |= gr_mat_mul_vector(tmp, G_top, H_top, e_last, ctx);
    status |= gr_mat_apply_Z(Za_elast, tmp, ctx);
    gr_mat_clear(tmp, ctx);
  }

  // calculate Z * a^T * e_{last} (bleed for c)
  gr_mat_init(ZaT_elast, n1, 1, ctx);
  {
    gr_mat_t tmp;
    gr_mat_init(tmp, n1, 1, ctx);
    status |= gr_mat_mul_vector(tmp, H_top, G_top, e_last, ctx);
    status |= gr_mat_apply_Z(ZaT_elast, tmp, ctx);
    gr_mat_clear(tmp, ctx);
  }

  // e0 of size n2 — appended to H_b
  // (Hb has n2 rows, so this vector must have n2 rows)
  gr_mat_init(e0_n2_Hb, n2, 1, ctx);
  status |= gr_mat_zero(e0_n2_Hb, ctx);
  status |= gr_one(gr_mat_entry_ptr(e0_n2_Hb, 0, 0, ctx), ctx);

  // e0 of size n2 — appended to G_c
  // (Gc has n2 rows, so this vector must have n2 rows)
  gr_mat_init(e0_n2_Gc, n2, 1, ctx);
  status |= gr_mat_zero(e0_n2_Gc, ctx);
  status |= gr_one(gr_mat_entry_ptr(e0_n2_Gc, 0, 0, ctx), ctx);

  if (status != GR_SUCCESS) {
    printf("[inv] FAILED: building correction vectors\n");
    gr_mat_clear(e_last, ctx);
    gr_mat_clear(Za_elast, ctx);
    gr_mat_clear(ZaT_elast, ctx);
    gr_mat_clear(e0_n2_Hb, ctx);
    gr_mat_clear(e0_n2_Gc, ctx);
    goto free_windows;
  }

  /* ASSEMBLING THE QUADRANTS: -------------------------------------------
   * We can now build our generators for these 4 quadrants.
   * Block a is just G_top, H_top
   * Block b combines G_top H_bottom with bleed correction
   * Block c combines G_bottom H_top with bleed correction
   * Block d combines the bottom generators and the top generators
   */

  // 2.2: b : Gb = [Gtop | Za_elast],  Hb = [Hbottom | e0]
  gr_mat_t G_b, H_b;
  slong rank_b = rank + 1;
  gr_mat_init(G_b, n1, rank_b, ctx);
  gr_mat_init(H_b, n2, rank_b, ctx);
  status |= gr_mat_concat_horizontal(G_b, G_top, Za_elast, ctx);
  status |= gr_mat_concat_horizontal(H_b, H_bottom, e0_n2_Hb, ctx);
  if (status != GR_SUCCESS)
    printf("[inv] FAILED: building G_b/H_b\n");
  else {
    printf("[inv] G_b:\n");
    gr_mat_print(G_b, ctx);
    printf("[inv] H_b:\n");
    gr_mat_print(H_b, ctx);
  }
  gr_mat_clear(Za_elast, ctx);
  gr_mat_clear(e0_n2_Hb, ctx);

  // 2.3: c : Gc = [Gbottom | e0],  Hc = [Htop | ZaT_elast]
  gr_mat_t G_c, H_c;
  slong rank_c = rank + 1;
  gr_mat_init(G_c, n2, rank_c, ctx);
  gr_mat_init(H_c, n1, rank_c, ctx);
  status |= gr_mat_concat_horizontal(G_c, G_bottom, e0_n2_Gc, ctx);
  status |= gr_mat_concat_horizontal(H_c, H_top, ZaT_elast, ctx);
  if (status != GR_SUCCESS)
    printf("[inv] FAILED: building G_c/H_c\n");
  else {
    printf("[inv] G_c:\n");
    gr_mat_print(G_c, ctx);
    printf("[inv] H_c:\n");
    gr_mat_print(H_c, ctx);
  }
  gr_mat_clear(ZaT_elast, ctx);
  gr_mat_clear(e0_n2_Gc, ctx);
  gr_mat_clear(e_last, ctx);

  // 2.4: d : Gd = [Gbottom | Gtop],  Hd = [Hbottom | Htop]
  gr_mat_t G_d, H_d;
  status |= gr_mat_addition_generateur(G_bottom, H_bottom, G_top, H_top, G_d, H_d, ctx);
  if (status != GR_SUCCESS) {
    printf("[inv] FAILED: building G_d/H_d\n");
    gr_mat_clear(G_d, ctx);
    gr_mat_clear(H_d, ctx);
    gr_mat_clear(G_c, ctx);
    gr_mat_clear(H_c, ctx);
    gr_mat_clear(G_b, ctx);
    gr_mat_clear(H_b, ctx);
    goto free_windows;
  }
  printf("[inv] G_d:\n");
  gr_mat_print(G_d, ctx);
  printf("[inv] H_d:\n");
  gr_mat_print(H_d, ctx);

  /* RECURSION #1 - INVERT BLOCK a -------------------------------------------
   * To apply the strassen's block inversion we need the inverse of the block a.
   * We call the fuction on the generators of the block a. The result is placed on G_e and H_e.
   */

  // 3: calculate rec e := a^{-1}
  gr_mat_t G_e, H_e;
  gr_mat_init(G_e, n1, rank, ctx);
  gr_mat_init(H_e, n1, rank, ctx);
  printf("[inv] --- recursive call for e=a^{-1} ---\n");
  status |= gr_toeplitz_inverse(G_e, H_e, G_top, H_top, ctx);
  if (status != GR_SUCCESS) {
    printf("[inv] FAILED: recursive inverse of a\n");
    goto free_e_top_bottom;
  }
  printf("[inv] G_e:\n");
  gr_mat_print(G_e, ctx);
  printf("[inv] H_e:\n");
  gr_mat_print(H_e, ctx);

  /* SHUR COMPLEMENT S=d−c*e*b -------------------------------------------
   * As we multiply these generators their ranks explose.
   * We call the compression algorithm after each arithmetics to keep
   * the algorithm efficient.
   *
   * We do 3 multiplications in sequence.
   */

  // 4: calculate generators for S := d - c*e*b

  // 4.1 ce = c * e
  gr_mat_t G_ce, H_ce;
  status |= gr_mat_mul_generator(G_ce, H_ce, G_c, H_c, G_e, H_e, ctx);
  status |= gr_mat_generator_compress(G_ce, H_ce, ctx);
  gr_mat_clear(G_c, ctx);
  gr_mat_clear(H_c, ctx);
  if (status != GR_SUCCESS) {
    printf("[inv] FAILED: computing ce = c*e\n");
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_e_top_bottom;
  }
  printf("[inv] G_ce (rank=%ld):\n", gr_mat_ncols(G_ce, ctx));
  gr_mat_print(G_ce, ctx);

  // 4.2 eb = e * b
  gr_mat_t G_eb, H_eb;
  status |= gr_mat_mul_generator(G_eb, H_eb, G_e, H_e, G_b, H_b, ctx);
  status |= gr_mat_generator_compress(G_eb, H_eb, ctx);
  gr_mat_clear(G_b, ctx);
  gr_mat_clear(H_b, ctx);
  if (status != GR_SUCCESS) {
    printf("[inv] FAILED: computing eb = e*b\n");
    gr_mat_clear(G_eb, ctx);
    gr_mat_clear(H_eb, ctx);
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_e_top_bottom;
  }
  printf("[inv] G_eb (rank=%ld):\n", gr_mat_ncols(G_eb, ctx));
  gr_mat_print(G_eb, ctx);

  // 4.3 ceb = ce * eb, then negate G in-place
  gr_mat_t G_ceb, H_ceb;
  status |= gr_mat_mul_generator(G_ceb, H_ceb, G_ce, H_ce, G_eb, H_eb, ctx);
  status |= gr_mat_neg(G_ceb, G_ceb, ctx); // Negate to get -ceb
  status |= gr_mat_generator_compress(G_ceb, H_ceb, ctx);
  if (status != GR_SUCCESS) {
    printf("[inv] FAILED: computing ceb = ce*eb\n");
    gr_mat_clear(G_ceb, ctx);
    gr_mat_clear(H_ceb, ctx);
    gr_mat_clear(G_eb, ctx);
    gr_mat_clear(H_eb, ctx);
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_e_top_bottom;
  }
  printf("[inv] G_ceb (rank=%ld):\n", gr_mat_ncols(G_ceb, ctx));
  gr_mat_print(G_ceb, ctx);

  // 4.4 S = d + (-ceb)
  gr_mat_t G_S, H_S;
  status |= gr_mat_addition_generateur(G_d, H_d, G_ceb, H_ceb, G_S, H_S, ctx);
  status |= gr_mat_generator_compress(G_S, H_S, ctx);
  gr_mat_clear(G_ceb, ctx);
  gr_mat_clear(H_ceb, ctx);
  gr_mat_clear(G_d, ctx);
  gr_mat_clear(H_d, ctx);
  if (status != GR_SUCCESS) {
    printf("[inv] FAILED: computing S = d - ceb\n");
    gr_mat_clear(G_S, ctx);
    gr_mat_clear(H_S, ctx);
    gr_mat_clear(G_eb, ctx);
    gr_mat_clear(H_eb, ctx);
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_e_top_bottom;
  }
  printf("[inv] G_S (rank=%ld):\n", gr_mat_ncols(G_S, ctx));
  gr_mat_print(G_S, ctx);
  printf("[inv] H_S:\n");
  gr_mat_print(H_S, ctx);

  /* RECURSION #2 - INVERT SHUR COMPLEMENT -------------------------------------------
   * Algoritm calls the inversion of the S.
   * We call the fuction on the generators of S. The result is placed on G_t and H_t.
   */

  // 5: calculate inv generators for t := S^{-1}
  gr_mat_t G_t, H_t;
  gr_mat_init(G_t, n2, gr_mat_ncols(G_S, ctx), ctx);
  gr_mat_init(H_t, n2, gr_mat_ncols(H_S, ctx), ctx);
  printf("[inv] --- recursive call for t=S^{-1} ---\n");
  status |= gr_toeplitz_inverse(G_t, H_t, G_S, H_S, ctx);
  gr_mat_clear(G_S, ctx);
  gr_mat_clear(H_S, ctx);
  if (status != GR_SUCCESS) {
    printf("[inv] FAILED: recursive inverse of S\n");
    gr_mat_clear(G_eb, ctx);
    gr_mat_clear(H_eb, ctx);
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_t;
  }
  printf("[inv] G_t:\n");
  gr_mat_print(G_t, ctx);
  printf("[inv] H_t:\n");
  gr_mat_print(H_t, ctx);

  /* FINAL STRASSEN ASSEMBLY -------------------------------------------
   * We have all our base components (e,eb,ce,t). We now assemble the four
   * quadrants of our final matrix A^(-1) [[x,y],[z,t]] where the algorithm
   * states x = e + ebtce, y = -ebt and z = -tce
   */

  // 6: calculate inv generators for A^(-1) [[x,y],[z,t]]
  // ebt = eb * t
  gr_mat_t G_ebt, H_ebt;
  status |= gr_mat_mul_generator(G_ebt, H_ebt, G_eb, H_eb, G_t, H_t, ctx);
  status |= gr_mat_generator_compress(G_ebt, H_ebt, ctx);
  gr_mat_clear(G_eb, ctx);
  gr_mat_clear(H_eb, ctx);
  if (status != GR_SUCCESS) {
    printf("[inv] FAILED: computing ebt = eb*t\n");
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_ebt;
  }

  // tce = t * ce
  gr_mat_t G_tce, H_tce;
  status |= gr_mat_mul_generator(G_tce, H_tce, G_t, H_t, G_ce, H_ce, ctx);
  status |= gr_mat_generator_compress(G_tce, H_tce, ctx);
  if (status != GR_SUCCESS) {
    printf("[inv] FAILED: computing tce = t*ce\n");
    gr_mat_clear(G_ce, ctx);
    gr_mat_clear(H_ce, ctx);
    goto free_tce;
  }

  // ebtce = ebt * ce
  gr_mat_t G_ebtce, H_ebtce;
  status |= gr_mat_mul_generator(G_ebtce, H_ebtce, G_ebt, H_ebt, G_ce, H_ce, ctx);
  status |= gr_mat_generator_compress(G_ebtce, H_ebtce, ctx);
  gr_mat_clear(G_ce, ctx);
  gr_mat_clear(H_ce, ctx);
  if (status != GR_SUCCESS) {
    printf("[inv] FAILED: computing ebtce = ebt*ce\n");
    gr_mat_clear(G_ebtce, ctx);
    gr_mat_clear(H_ebtce, ctx);
    goto free_tce;
  }

  // x = e + ebtce
  gr_mat_t G_x, H_x;
  status |= gr_mat_addition_generateur(G_e, H_e, G_ebtce, H_ebtce, G_x, H_x, ctx);
  status |= gr_mat_generator_compress(G_x, H_x, ctx);
  gr_mat_clear(G_ebtce, ctx);
  gr_mat_clear(H_ebtce, ctx);
  if (status != GR_SUCCESS) {
    printf("[inv] FAILED: computing x = e + ebtce\n");
    goto free_x;
  }
  printf("[inv] G_x:\n");
  gr_mat_print(G_x, ctx);

  // y = -ebt (negate G in-place)
  status |= gr_mat_neg(G_ebt, G_ebt, ctx);
  if (status != GR_SUCCESS) {
    printf("[inv] FAILED: negate G_ebt\n");
    goto free_x;
  }

  // z = -tce (negate G in-place)
  status |= gr_mat_neg(G_tce, G_tce, ctx);
  if (status != GR_SUCCESS) {
    printf("[inv] FAILED: negate G_tce\n");
    goto free_x;
  }
  
  /* PACKING THE OUTPUT -------------------------------------------
   * We have now our generators for x, y, z and t.
   * Since they represent blocks in a larger nxn matrix, their global generators
   * are just these individual generators placed in thier correct row offsets.
   * 
   * We create an empty n x total_rank matrices of G_D and H_D. We loop through 
   * x, y, z and t into their correct rows (top half for n1 or bottom half for n2) 
   * Then we compress one last time.
   */

  // pack [x y; z t] into G_D / H_D
  {
    slong rx = gr_mat_ncols(G_x, ctx);
    slong ry = gr_mat_ncols(G_ebt, ctx);
    slong rz = gr_mat_ncols(G_tce, ctx);
    slong rt = gr_mat_ncols(G_t, ctx);
    slong total_rank = rx + ry + rz + rt;
    printf("[inv] packing: rx=%ld ry=%ld rz=%ld rt=%ld\n", rx, ry, rz, rt);

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
    if (status != GR_SUCCESS) printf("[inv] FAILED: packing G_D/H_D\n");
  }

  status |= gr_mat_generator_compress(G_D, H_D, ctx);
  if (status != GR_SUCCESS)
    printf("[inv] FAILED: final compress G_D/H_D\n");
  else {
    printf("[inv] G_D:\n");
    gr_mat_print(G_D, ctx);
    printf("[inv] H_D:\n");
    gr_mat_print(H_D, ctx);
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