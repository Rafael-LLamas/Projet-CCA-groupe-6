#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"
#include "flint/ulong_extras.h"
#include "gr_struct_mat.h"
#include "util.h"

int test_addition_generateurs_toeplitz(int nb_iter)
{
  int i = 0;
  int status = GR_SUCCESS;
  flint_rand_t state;
  flint_rand_init(state);
  while (i < nb_iter && status == GR_SUCCESS)
  {
    gr_struct_mat_t A, B;
    gr_mat_t A_DENSE, B_DENSE, C_DENSE, CHECK_DENSE;
    gr_ctx_t ctx;

    flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
    gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

    slong n = n_randint(state, 100);
    slong m = n_randint(state, 100);

    gr_struct_mat_init(A, n, m, T_TOEPLITZ, DISP_PLUS, ctx);
    gr_struct_mat_init(B, n, m, T_TOEPLITZ, DISP_PLUS, ctx);
    gr_mat_init(A_DENSE, n, m, ctx);
    gr_mat_init(B_DENSE, n, m, ctx);
    gr_mat_init(C_DENSE, n, m, ctx);
    gr_mat_init(CHECK_DENSE, n, m, ctx);

    status |= gr_struct_mat_random(A, state, ctx);
    status |= gr_struct_mat_random(B, state, ctx);

    status |= gr_struct_mat_reconstruct(A_DENSE, A, ctx);
    status |= gr_struct_mat_reconstruct(B_DENSE, B, ctx);

    status |= gr_mat_add(C_DENSE, A_DENSE, B_DENSE, ctx);

    status |= gr_struct_mat_add(A, B, ctx); // A = A + B

    status |= gr_struct_mat_reconstruct(CHECK_DENSE, A, ctx);

    if (gr_mat_equal(C_DENSE, CHECK_DENSE, ctx) != T_TRUE) status = GR_TEST_FAIL;

    gr_struct_mat_clear(A, ctx);
    gr_struct_mat_clear(B, ctx);
    gr_mat_clear(A_DENSE, ctx);
    gr_mat_clear(B_DENSE, ctx);
    gr_mat_clear(C_DENSE, ctx);
    gr_mat_clear(CHECK_DENSE, ctx);
    gr_ctx_clear(ctx);
    i++;
  }

  flint_rand_clear(state);
  return status;
}

// int test_addition_generateurs_quasi_toeplitz()
// {
//   int i = 0;
//   int status = GR_SUCCESS;
//   flint_rand_t state;
//   flint_rand_init(state);
//   while (i < 50 && status == GR_SUCCESS)
//   {
//     gr_mat_t A, B, C;
//     gr_ctx_t ctx;
//     flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
//     gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));
//     slong n = n_randint(state, 100);
//     slong m = n_randint(state, 100);
//     slong r = n_randint(state, 5);

//     gr_mat_init(A, n, m, ctx);
//     gr_mat_init(B, n, m, ctx);
//     gr_mat_init(C, n, m, ctx);

//     status = gr_mat_quasi_toeplitz_rank(A, r, state, ctx);

//     status = gr_mat_quasi_toeplitz_rank(B, r, state, ctx);

//     status = gr_mat_add(C, A, B, ctx);

//     gr_mat_t A_G, A_H, B_G, B_H, D_G, D_H, D;
//     gr_mat_init(D, n, m, ctx);

//     status = gr_mat_G_H(A_G, A_H, A, DISP_PLUS, ctx);

//     status = gr_mat_G_H(B_G, B_H, B, DISP_PLUS, ctx);

//     status = gr_mat_addition_generateur(A_G, A_H, B_G, B_H, D_G, D_H, ctx);

//     status = gr_mat_reconstruct_A(D, D_G, D_H, DISP_PLUS, ctx);

//     if (gr_mat_equal(C, D, ctx) != T_TRUE) i = 9998;

//     gr_mat_clear(A, ctx);
//     gr_mat_clear(B, ctx);
//     gr_mat_clear(C, ctx);
//     gr_mat_clear(A_G, ctx);
//     gr_mat_clear(A_H, ctx);
//     gr_mat_clear(B_G, ctx);
//     gr_mat_clear(B_H, ctx);
//     gr_mat_clear(D_G, ctx);
//     gr_mat_clear(D_H, ctx);
//     gr_mat_clear(D, ctx);
//     gr_ctx_clear(ctx);
//     i++;
//   }
//   flint_rand_clear(state);
//   return status;
// }

void usage(char *argv[])
{
  fprintf(stderr, "Usage: %s <test_name>\n", argv[0]);
  fprintf(stderr, "Available tests:\n");
  fprintf(stderr, "  - addition_generators_toeplitz\n");
  fprintf(stderr, "  - addition_generators_quasi_toeplitz\n");
}
int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    usage(argv);
    return GR_UNABLE;
  }

  fprintf(stderr, "=> Start test \"%s\"\n", argv[1]);
  int ok = GR_SUCCESS;

  if (strcmp("addition_generators_toeplitz", argv[1]) == 0) { ok = test_addition_generateurs_toeplitz(20); }
  // else if (strcmp("addition_generators_quasi_toeplitz", argv[1]) == 0) { ok = test_addition_generateurs_quasi_toeplitz(); }
  else
  {
    fprintf(stderr, "status: test \"%s\" not found!\n", argv[1]);
    exit(EXIT_FAILURE);
  }

  if (ok == GR_SUCCESS)
  {
    fprintf(stderr, "Test \"%s\" finished: SUCCESS\n", argv[1]);
    return EXIT_SUCCESS;
  }
  else
  {
    fprintf(stderr, "Test \"%s\" finished: FAILURE\n", argv[1]);
    return EXIT_FAILURE;
  }
}