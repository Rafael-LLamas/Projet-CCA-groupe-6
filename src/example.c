#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"
#include "gr_struct_mat.h"
#include <stdio.h>
#include <stdlib.h>

#define N 4

void print_mat(const char *label, gr_mat_t mat, gr_ctx_t ctx)
{
  printf("%s\n", label);
  gr_mat_print(mat, ctx);
  printf("\n");
}


int example_toeplitz(gr_ctx_t ctx)
{
  int status = GR_SUCCESS;

  printf("========================================\n");
  printf("TOEPLITZ EXAMPLE (DISP_PLUS)\n");
  printf("========================================\n\n");

  // build the dense Toeplitz matrix
  gr_mat_t src;
  gr_mat_init(src, N, N, ctx);
  for (slong i = 0; i < N; i++)
    for (slong j = 0; j < N; j++) status |= gr_set_si(gr_mat_entry_ptr(src, i, j, ctx), i - j, ctx);
  print_mat("Source Toeplitz matrix:", src, ctx);

  gr_struct_mat_t smat;
  status |= gr_struct_mat_init_set(smat, src, T_TOEPLITZ, DISP_PLUS, ctx);

  printf("Structured representation:\n");
  gr_struct_mat_print(smat, ctx);

  // reconstruct and compare
  gr_mat_t rec;
  gr_mat_init(rec, N, N, ctx);
  status |= gr_struct_mat_reconstruct(rec, smat, ctx);
  print_mat("Reconstructed matrix:", rec, ctx);


  gr_struct_mat_clear(smat, ctx);
  gr_mat_clear(rec, ctx);
  gr_mat_clear(src, ctx);

  printf("Toeplitz status: %s\n\n", status == GR_SUCCESS ? "GR_SUCCESS" : "FAIL");
  return status;
}

int example_hankel(gr_ctx_t ctx)
{
  int status = GR_SUCCESS;

  printf("========================================\n");
  printf("HANKEL EXAMPLE (DISP_MINUS)\n");
  printf("========================================\n\n");

  // build the dense Hankel matrix
  gr_mat_t src;
  gr_mat_init(src, N, N, ctx);
  for (slong i = 0; i < N; i++)
    for (slong j = 0; j < N; j++) status |= gr_set_si(gr_mat_entry_ptr(src, i, j, ctx), i + j + 1, ctx);
  print_mat("Source Hankel matrix:", src, ctx);

  gr_struct_mat_t smat;
  status |= gr_struct_mat_init_set(smat, src, T_HANKEL, DISP_MINUS, ctx);
  printf("Structured representation:\n");
  gr_struct_mat_print(smat, ctx);

  // reconstruct and compare
  gr_mat_t rec;
  gr_mat_init(rec, N, N, ctx);
  status |= gr_struct_mat_reconstruct(rec, smat, ctx);
  print_mat("Reconstructed matrix:", rec, ctx);

  gr_struct_mat_clear(smat, ctx);
  gr_mat_clear(rec, ctx);
  gr_mat_clear(src, ctx);

  printf("Hankel status: %s\n\n", status == GR_SUCCESS ? "OK" : "FAIL");
  return status;
}

int main(void)
{
  int status = GR_SUCCESS;

  gr_ctx_t ctx;
  gr_ctx_init_fmpq(ctx);

  status |= example_toeplitz(ctx);
  status |= example_hankel(ctx);

  gr_ctx_clear(ctx);
  return status;
}