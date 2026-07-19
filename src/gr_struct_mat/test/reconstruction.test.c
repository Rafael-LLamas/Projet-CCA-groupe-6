#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "gr_struct_mat.h"
#include "utility.h"

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/ulong_extras.h"

// Checks D(A) == G * H^T
int verify_reconstruction(gr_mat_t A, gr_struct_mat_t S, gr_ctx_t ctx)
{
  int status = GR_SUCCESS;

  slong n = gr_mat_nrows(A, ctx);
  slong m = gr_mat_ncols(A, ctx);
  slong rank = gr_mat_ncols(S->G, ctx);

  gr_mat_t D, HT, GHT;
  gr_mat_init(D, n, m, ctx);
  gr_mat_init(HT, rank, m, ctx);
  gr_mat_init(GHT, n, m, ctx);

  status |= _gr_mat_displacement(D, A, S->struct_t, S->disp_t, ctx);

  status |= gr_mat_transpose(HT, S->H, ctx);
  status |= gr_mat_mul(GHT, S->G, HT, ctx);

  if (gr_mat_equal(D, GHT, ctx) == T_FALSE)
  {
    flint_printf("D(A) != G * H^T\n");
    status = GR_TEST_FAIL;
  }

  gr_mat_clear(D, ctx);
  gr_mat_clear(HT, ctx);
  gr_mat_clear(GHT, ctx);

  return status;
}


int test_toeplitz_reconstruction(int nb_iter)
{
  flint_rand_t state;
  flint_rand_init(state);
  int status = GR_SUCCESS;

  for (int i = 0; i < nb_iter; i++)
  {
    ulong prime = n_randtest_prime(state, 0);
    gr_ctx_t ctx;
    gr_ctx_init_nmod(ctx, prime);

    slong n = n_randint(state, 100) + 10;
    slong m = n_randint(state, 100) + 10;

    gr_struct_mat_t S;
    gr_mat_t A;

    gr_struct_mat_init(S, n, m, T_TOEPLITZ, DISP_PLUS, ctx);
    status |= gr_struct_mat_random(S, state, ctx);
    if (status != GR_SUCCESS) break;

    gr_mat_init(A, n, m, ctx);
    status |= gr_struct_mat_reconstruct(A, S, ctx);
    if (status != GR_SUCCESS) break;

    status |= verify_reconstruction(A, S, ctx);
    if (status != GR_SUCCESS)
    {
      flint_printf("Failed on pure Toeplitz iteration %d\n", i);
      break;
    }

    gr_struct_mat_clear(S, ctx);
    gr_mat_clear(A, ctx);
    gr_ctx_clear(ctx);
  }

  flint_rand_clear(state);
  return status;
}

int test_quasi_toeplitz_reconstruction(int nb_iter)
{
  flint_rand_t state;
  flint_rand_init(state);
  int status = GR_SUCCESS;

  for (int i = 0; i < nb_iter; i++)
  {
    ulong prime = n_randtest_prime(state, 0);
    gr_ctx_t ctx;
    gr_ctx_init_nmod(ctx, prime);

    slong n = n_randint(state, 100) + 10;
    slong m = n_randint(state, 100) + 10;

    gr_struct_mat_t S1, S2;
    gr_mat_t A;

    gr_struct_mat_init(S1, n, m, T_TOEPLITZ, DISP_PLUS, ctx);
    gr_struct_mat_init(S2, n, m, T_TOEPLITZ, DISP_PLUS, ctx);

    status |= gr_struct_mat_random(S1, state, ctx);
    status |= gr_struct_mat_random(S2, state, ctx);

    // S1 is quasi-toeplitz
    status |= gr_struct_mat_add(S1, S2, ctx);
    if (status != GR_SUCCESS) break;

    gr_mat_init(A, n, m, ctx);
    status |= gr_struct_mat_reconstruct(A, S1, ctx);
    if (status != GR_SUCCESS) break;

    status |= verify_reconstruction(A, S1, ctx);
    if (status != GR_SUCCESS)
    {
      flint_printf("Failed on Quasi-Toeplitz iteration %d\n", i);
      break;
    }

    gr_struct_mat_clear(S1, ctx);
    gr_struct_mat_clear(S2, ctx);
    gr_mat_clear(A, ctx);
    gr_ctx_clear(ctx);
  }

  flint_rand_clear(state);
  return status;
}

void usage(char *argv[]) { fprintf(stderr, "Usage: %s <test_name>\n", argv[0]); }

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    usage(argv);
    return GR_UNABLE;
  }

  fprintf(stderr, "=> Start test \"%s\"\n", argv[1]);
  int ok = GR_SUCCESS;

  if (strcmp("toeplitz_reconstruction", argv[1]) == 0) { ok = test_toeplitz_reconstruction(30); }
  else if (strcmp("quasi_toeplitz_reconstruction", argv[1]) == 0) { ok = test_quasi_toeplitz_reconstruction(30); }
  else
  {
    fprintf(stderr, "status: test \"%s\" not found!\n", argv[1]);
    exit(EXIT_FAILURE);
  }

  if (!ok)
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