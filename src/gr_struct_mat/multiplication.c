#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "utility/utility.h"
#include <stdlib.h>
#include <time.h>

// Multiply a vector with a matrix represented by its generators
int _gr_mat_mul_vector_toeplitz(gr_mat_t Res, gr_mat_t G, gr_mat_t H, gr_mat_t V, gr_ctx_t ctx)
{
  int status = GR_SUCCESS;
  slong n = gr_mat_nrows(G, ctx), m = gr_mat_nrows(H, ctx);
  slong alpha = gr_mat_ncols(G, ctx), m_cols = gr_mat_ncols(V, ctx);

  gr_poly_t pg, ph, pv, p_tmp, p_tmp_shift, p_tmp2;

  status |= gr_mat_zero(Res, ctx);

  for (slong j = 0; j < m_cols; j++)
  {
    gr_poly_init(pv, ctx);
    gr_poly_fit_length(pv, m, ctx);
    for (slong r = 0; r < m; r++) { status |= gr_set(gr_poly_coeff_ptr(pv, r, ctx), gr_mat_entry_srcptr(V, r, j, ctx), ctx); }
    _gr_poly_set_length(pv, m, ctx);

    for (slong k = 0; k < alpha; k++)
    {
      gr_poly_init(pg, ctx);
      gr_poly_init(ph, ctx);
      gr_poly_init(p_tmp, ctx);
      gr_poly_init(p_tmp_shift, ctx);
      gr_poly_init(p_tmp2, ctx);
      gr_poly_fit_length(pg, n, ctx);
      for (slong r = 0; r < n; r++) { status |= gr_set(gr_poly_coeff_ptr(pg, r, ctx), gr_mat_entry_srcptr(G, r, k, ctx), ctx); }
      _gr_poly_set_length(pg, n, ctx);

      gr_poly_fit_length(ph, m, ctx);
      for (slong r = 0; r < m; r++) { status |= gr_set(gr_poly_coeff_ptr(ph, r, ctx), gr_mat_entry_srcptr(H, r, k, ctx), ctx); }
      _gr_poly_set_length(ph, m, ctx);

      status |= gr_poly_reverse(ph, ph, m, ctx);
      status |= gr_poly_mul(p_tmp, ph, pv, ctx);
      status |= gr_poly_shift_right(p_tmp_shift, p_tmp, m - 1, ctx);
      status |= gr_poly_mul(p_tmp2, pg, p_tmp_shift, ctx);

      if (status != GR_SUCCESS)
      {
        gr_poly_clear(pg, ctx);
        gr_poly_clear(ph, ctx);
        gr_poly_clear(pv, ctx);
        gr_poly_clear(p_tmp, ctx);
        gr_poly_clear(p_tmp_shift, ctx);
        gr_poly_clear(p_tmp2, ctx);
        return status;
      }

      if (gr_poly_length(p_tmp2, ctx) > n) _gr_poly_set_length(p_tmp2, n, ctx);

      for (slong i = 0; i < gr_poly_length(p_tmp2, ctx); i++)
      {
        status |= gr_add(gr_mat_entry_ptr(Res, i, j, ctx), gr_mat_entry_srcptr(Res, i, j, ctx), gr_poly_coeff_srcptr(p_tmp2, i, ctx), ctx);
      }
      gr_poly_clear(pg, ctx);
      gr_poly_clear(ph, ctx);
      gr_poly_clear(p_tmp_shift, ctx);
      gr_poly_clear(p_tmp, ctx);
      gr_poly_clear(p_tmp2, ctx);
    }
    gr_poly_clear(pv, ctx);
  }

  return status;
}

// Multiply two matrices represented by generators to a third destination represented by generators
int _gr_mat_mul_generator_toeplitz(gr_mat_t G_c, gr_mat_t H_c, gr_mat_t T, gr_mat_t U, gr_mat_t G, gr_mat_t H, gr_ctx_t ctx)
{
  slong n = gr_mat_nrows(T, ctx); // Lignes de A
  slong m = gr_mat_nrows(G, ctx); // Colonnes de A
  slong k = gr_mat_nrows(H, ctx); // colonnes de B

  gr_mat_t W, V, a, b, Tmp, LastCol_m, LastCol_k;
  int status = GR_SUCCESS;

  // Initialisation avec les dimensions respectives
  gr_mat_init(W, n, gr_mat_ncols(G, ctx), ctx);
  gr_mat_init(V, k, gr_mat_ncols(U, ctx), ctx);
  gr_mat_init(a, n, 1, ctx);
  gr_mat_init(b, k, 1, ctx);
  gr_mat_init(Tmp, m, gr_mat_ncols(G, ctx), ctx);

  // 1. W = Z * A * Z^T * G
  status |= _gr_mat_apply_Zt(Tmp, G, ctx);
  status |= _gr_mat_mul_vector_toeplitz(W, T, U, Tmp, ctx);
  status |= _gr_mat_apply_Z(W, W, ctx);


  // 2. V = B^T * U (Le générateur de B^T est {H, G})
  status |= _gr_mat_mul_vector_toeplitz(V, H, G, U, ctx);

  // 3. Correction colonnes
  // Vecteur e_{m-1} pour la dimension de A
  // LastCol_m doit être e_{nrows(U) - 1}
  slong m_a = gr_mat_nrows(U, ctx);
  gr_mat_init(LastCol_m, m_a, 1, ctx);

  status |= gr_mat_zero(LastCol_m, ctx);
  status |= gr_one(gr_mat_entry_ptr(LastCol_m, m_a - 1, 0, ctx), ctx);
  status |= _gr_mat_mul_vector_toeplitz(a, T, U, LastCol_m, ctx);
  status |= _gr_mat_apply_Z(a, a, ctx);

  // Vecteur e_{m-1} pour la dimension de B^T
  // B^T a pour générateurs {H, G}, donc m = nrows(G)
  slong m_bt = gr_mat_nrows(G, ctx);
  gr_mat_init(LastCol_k, m_bt, 1, ctx);

  status |= gr_mat_zero(LastCol_k, ctx);
  status |= gr_one(gr_mat_entry_ptr(LastCol_k, m_bt - 1, 0, ctx), ctx);
  status |= _gr_mat_mul_vector_toeplitz(b, H, G, LastCol_k, ctx);
  status |= _gr_mat_apply_Z(b, b, ctx);

  gr_mat_clear(Tmp, ctx);
  gr_mat_clear(LastCol_m, ctx);
  gr_mat_clear(LastCol_k, ctx);

  // 4. Concaténation (G_c: n x (ra+rb+1), H_c: k x (ra+rb+1))
  gr_mat_t G_temp, H_temp;
  gr_mat_init(G_c, n, gr_mat_ncols(T, ctx) + gr_mat_ncols(W, ctx) + 1, ctx);
  gr_mat_init(G_temp, n, gr_mat_ncols(T, ctx) + gr_mat_ncols(W, ctx), ctx);
  gr_mat_init(H_c, k, gr_mat_ncols(V, ctx) + gr_mat_ncols(H, ctx) + 1, ctx);
  gr_mat_init(H_temp, k, gr_mat_ncols(V, ctx) + gr_mat_ncols(H, ctx), ctx);

  status |= gr_mat_concat_horizontal(G_temp, T, W, ctx);
  status |= gr_mat_concat_horizontal(G_c, G_temp, a, ctx);
  status |= gr_mat_neg(b, b, ctx);
  status |= gr_mat_concat_horizontal(H_temp, V, H, ctx);
  status |= gr_mat_concat_horizontal(H_c, H_temp, b, ctx);

  gr_mat_clear(W, ctx);
  gr_mat_clear(V, ctx);
  gr_mat_clear(a, ctx);
  gr_mat_clear(b, ctx);
  gr_mat_clear(G_temp, ctx);
  gr_mat_clear(H_temp, ctx);

  return status;
}


// (==============================================================)


// Does mat1 * mat2 = dest, dest should be initialized.
int gr_struct_mat_mul(gr_struct_mat_t dest, gr_struct_mat_t mat1, gr_struct_mat_t mat2, gr_ctx_t ctx)
{
  if (gr_struct_mat_ncols(mat1, ctx) != gr_struct_mat_ncols(dest, ctx) || gr_struct_mat_nrows(mat2, ctx) != gr_struct_mat_nrows(dest, ctx))
    return GR_DOMAIN;

  int status = GR_SUCCESS;
  structure_type_t struct_t_var = mat1->struct_t;

  if (struct_t_var == T_TOEPLITZ)
    status |= _gr_mat_mul_generator_toeplitz(dest->G, dest->H, mat1->G, mat1->H, mat2->G, mat2->H, ctx);
  else
    status = GR_UNABLE;

  return status;
}


// Does dest * mat and replaces on dest
int gr_struct_mat_mul_f(gr_struct_mat_t dest, gr_struct_mat_t mat, gr_ctx_t ctx)
{
  if (gr_struct_mat_ncols(dest, ctx) != gr_struct_mat_nrows(mat, ctx)) return GR_DOMAIN;

  int status = GR_SUCCESS;
  structure_type_t struct_t_var = dest->struct_t;

  gr_struct_mat_t s_mat;
  gr_struct_mat_init(s_mat, gr_struct_mat_nrows(dest, ctx), gr_struct_mat_ncols(mat, ctx), struct_t_var, dest->disp_t, ctx);

  // TODO: Hankel Multiplication
  if (struct_t_var == T_TOEPLITZ) { status |= _gr_mat_mul_generator_toeplitz(s_mat->G, s_mat->H, dest->G, dest->H, mat->G, mat->H, ctx); }
  else
    status = GR_UNABLE;

  if (status == GR_SUCCESS)
  {
    gr_mat_swap(s_mat->G, dest->G, ctx);
    gr_mat_swap(s_mat->H, dest->H, ctx);
  }

  gr_struct_mat_clear(s_mat, ctx);
  return status;
}
