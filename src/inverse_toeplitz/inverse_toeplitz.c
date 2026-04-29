#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"

#include "addition.h"
#include "compression.h"
#include "flint/gr_types.h"
#include "inverse_toeplitz_aux.h"
#include "multiplication.h"

int gr_mat_inverse_toeplitz(gr_mat_t G_D, gr_mat_t H_D, gr_mat_t G_A, gr_mat_t H_A, gr_ctx_t ctx) {
  int status = GR_SUCCESS;
  slong n = gr_mat_nrows(G_A, ctx);
  slong rank = gr_mat_ncols(G_A, ctx);

  /* BASE CASE: -------------------------------------------
   * If a matrix is 1x1, its shift operator Z is 0.
   * Therefore its displacement DISP_PLUS and DISP_MINS is the same.
   *
   * Once recieved, we simply invert the value inside rthe G_A and
   * H_A and place it to G_D and H_D.
   */
  if (n == 1) {
    gr_ptr val = gr_heap_init(ctx);
    gr_ptr temp = gr_heap_init(ctx);

    status |= gr_zero(val, ctx);
    for (slong i = 0; i < rank; i++) {
      status |= gr_mul(temp, gr_mat_entry_srcptr(G_A, 0, i, ctx), gr_mat_entry_srcptr(H_A, 0, i, ctx), ctx);
      status |= gr_add(val, val, temp, ctx);
    }

    status = gr_inv(val, val, ctx);

    if (status == GR_SUCCESS) {
      gr_mat_clear(G_D, ctx);
      gr_mat_clear(H_D, ctx);
      gr_mat_init(G_D, 1, 1, ctx);
      gr_mat_init(H_D, 1, 1, ctx);

      status |= gr_set(gr_mat_entry_ptr(G_D, 0, 0, ctx), val, ctx);
      status |= gr_one(gr_mat_entry_ptr(H_D, 0, 0, ctx), ctx);
    }

    gr_heap_clear(val, ctx);
    gr_heap_clear(temp, ctx);
    return status;
  }

  /* GET BLOCKS a b c d ------------------------------------------- */
  gr_mat_t G_a, H_a, G_b, H_b, G_c, H_c, G_d, H_d;
  status |= gr_mat_split_quadrants(G_a, H_a, G_b, H_b, G_c, H_c, G_d, H_d, G_A, H_A, ctx);
  if (status != GR_SUCCESS) return status;

  /* RECURSION #1 - INVERT BLOCK a -------------------------------------------
   * To apply the strassen's block inversion we need the inverse of the block a.
   * We call the fuction on the generators of the block a. The result is placed on G_e and H_e.
   */
  gr_mat_t G_e, H_e;
  gr_mat_init(G_e, gr_mat_nrows(G_a, ctx), rank, ctx);
  gr_mat_init(H_e, gr_mat_nrows(H_a, ctx), rank, ctx);
  status |= gr_mat_inverse_toeplitz(G_e, H_e, G_a, H_a, ctx);

  /* SHUR COMPLEMENT S=d−c*e*b -------------------------------------------
   * As we multiply these generators their ranks explose.
   * We call the compression algorithm after each arithmetics to keep
   * the algorithm efficient.
   *
   * We do 3 multiplications in sequence.
   */
  gr_mat_t G_ce, H_ce, G_eb, H_eb, G_ceb, H_ceb, G_S, H_S;

  // ce = c * e
  status |= gr_mat_mul_generator(G_ce, H_ce, G_c, H_c, G_e, H_e, ctx);
  status |= gr_mat_generator_compress(G_ce, H_ce, ctx);

  // eb = e * b
  status |= gr_mat_mul_generator(G_eb, H_eb, G_e, H_e, G_b, H_b, ctx);
  status |= gr_mat_generator_compress(G_eb, H_eb, ctx);

  // ceb = c * eb
  status |= gr_mat_mul_generator(G_ceb, H_ceb, G_c, H_c, G_eb, H_eb, ctx);
  status |= gr_mat_neg(G_ceb, G_ceb, ctx); // -ceb
  status |= gr_mat_generator_compress(G_ceb, H_ceb, ctx);

  // S = d + (-ceb)
  status |= gr_mat_addition_generateur(G_d, H_d, G_ceb, H_ceb, G_S, H_S, ctx);
  status |= gr_mat_generator_compress(G_S, H_S, ctx);

  /* RECURSION #2 - INVERT SHUR COMPLEMENT -------------------------------------------
   * Algoritm calls the inversion of the S.
   * We call the fuction on the generators of S. The result is placed on G_t and H_t.
   */
  gr_mat_t G_t, H_t;
  gr_mat_init(G_t, gr_mat_nrows(G_S, ctx), gr_mat_ncols(G_S, ctx), ctx);
  gr_mat_init(H_t, gr_mat_nrows(H_S, ctx), gr_mat_ncols(H_S, ctx), ctx);
  status |= gr_mat_inverse_toeplitz(G_t, H_t, G_S, H_S, ctx);

  /* FINAL STRASSEN ASSEMBLY -------------------------------------------
   * We have all our base components (e,eb,ce,t). We now assemble the four
   * quadrants of our final matrix A^(-1) [[x,y],[z,t]] where the algorithm
   * states x = e + ebtce, y = -ebt and z = -tce
   */
  gr_mat_t G_ebt, H_ebt, G_tce, H_tce, G_ebtce, H_ebtce, G_x, H_x;

  // ebt = eb * t
  status |= gr_mat_mul_generator(G_ebt, H_ebt, G_eb, H_eb, G_t, H_t, ctx);
  status |= gr_mat_generator_compress(G_ebt, H_ebt, ctx);

  // tce = t * ce
  status |= gr_mat_mul_generator(G_tce, H_tce, G_t, H_t, G_ce, H_ce, ctx);
  status |= gr_mat_generator_compress(G_tce, H_tce, ctx);

  // ebtce = ebt * ce
  status |= gr_mat_mul_generator(G_ebtce, H_ebtce, G_ebt, H_ebt, G_ce, H_ce, ctx);
  status |= gr_mat_generator_compress(G_ebtce, H_ebtce, ctx);

  // x = e + ebtce
  status |= gr_mat_addition_generateur(G_e, H_e, G_ebtce, H_ebtce, G_x, H_x, ctx);
  status |= gr_mat_generator_compress(G_x, H_x, ctx);

  // y = -ebt
  status |= gr_mat_neg(G_ebt, G_ebt, ctx);

  // z = -tce
  status |= gr_mat_neg(G_tce, G_tce, ctx);

  /* PACK BLOCKS x y z t ------------------------------------------- */
  status |= gr_mat_pack_quadrants(G_D, H_D, G_x, H_x, G_ebt, H_ebt, G_tce, H_tce, G_t, H_t, ctx);

  gr_mat_clear(G_a, ctx);
  gr_mat_clear(H_a, ctx);
  gr_mat_clear(G_b, ctx);
  gr_mat_clear(H_b, ctx);
  gr_mat_clear(G_c, ctx);
  gr_mat_clear(H_c, ctx);
  gr_mat_clear(G_d, ctx);
  gr_mat_clear(H_d, ctx);

  gr_mat_clear(G_e, ctx);
  gr_mat_clear(H_e, ctx);
  gr_mat_clear(G_ce, ctx);
  gr_mat_clear(H_ce, ctx);
  gr_mat_clear(G_eb, ctx);
  gr_mat_clear(H_eb, ctx);
  gr_mat_clear(G_ceb, ctx);
  gr_mat_clear(H_ceb, ctx);
  gr_mat_clear(G_S, ctx);
  gr_mat_clear(H_S, ctx);
  gr_mat_clear(G_t, ctx);
  gr_mat_clear(H_t, ctx);

  gr_mat_clear(G_ebtce, ctx);
  gr_mat_clear(H_ebtce, ctx);
  gr_mat_clear(G_x, ctx);
  gr_mat_clear(H_x, ctx);
  gr_mat_clear(G_ebt, ctx);
  gr_mat_clear(H_ebt, ctx);
  gr_mat_clear(G_tce, ctx);
  gr_mat_clear(H_tce, ctx);

  return status;
}