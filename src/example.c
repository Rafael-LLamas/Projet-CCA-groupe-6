#include <stdio.h>
#include <stdlib.h>

#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"
#include "gr_struct_mat.h"

int main(void)
{
  int res = GR_SUCCESS;
  printf("Testing executable\n");

  gr_ctx_t ctx;
  gr_ctx_init_fmpz(ctx);

  // initialize the toeplitz matrix
  gr_mat_t testing_mat_source;
  gr_mat_init(testing_mat_source, 4, 4, ctx);

  for (slong i = 0; i < 4; i++)
  {
    for (slong j = 0; j < 4; j++)
    {
      slong toeplitz_val = i - j;
      res |= gr_set_si(gr_mat_entry_ptr(testing_mat_source, i, j, ctx), toeplitz_val, ctx);
    }
  }

  printf("\n\nSource Matrix:\n");
  gr_mat_print(testing_mat_source, ctx);
  printf("\n");

  // initialize the structured matrix structure
  gr_struct_mat_t testing_mat;

  res |= gr_struct_mat_init_set(testing_mat, testing_mat_source, T_TOEPLITZ, DISP_PLUS, ctx);
  // printf("\nStructured Matrix Type:\n");
  printf("\n");
  gr_struct_mat_print(testing_mat, ctx);

  gr_mat_t reconstructed_dense_mat;
  gr_mat_init(reconstructed_dense_mat, 0, 0, ctx);

  res |= gr_struct_mat_reconstruct(reconstructed_dense_mat, testing_mat, ctx);
  printf("\nMatrix reconstructed from structured matrix type:\n");
  gr_mat_print(reconstructed_dense_mat, ctx);
  printf("\n");

  gr_struct_mat_clear(testing_mat, ctx);
  gr_mat_clear(reconstructed_dense_mat, ctx);
  gr_mat_clear(testing_mat_source, ctx);
  gr_ctx_clear(ctx);
  return res;
}