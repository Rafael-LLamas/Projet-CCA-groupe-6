#include <stdio.h>
#include <stdlib.h>

#include "gr_struct_mat.h"

int main(void)
{
  printf("Testing executable\n");

  gr_ctx_t ctx;
  gr_ctx_init_fmpz(ctx);

  gr_struct_mat_t testing_mat;

  gr_struct_mat_init(testing_mat, 5, 5, T_TOEPLITZ, DISP_PLUS, ctx);
  gr_struct_mat_print(testing_mat, ctx);

  gr_struct_mat_clear(testing_mat, ctx);
  gr_ctx_clear(ctx);

  return EXIT_SUCCESS;
}