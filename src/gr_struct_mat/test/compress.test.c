#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "gr_struct_mat.h"

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/ulong_extras.h"

int verify_compression(gr_mat_t ref, gr_struct_mat_t S, slong max_expected_rank, gr_ctx_t ctx)
{
  int status = GR_SUCCESS;
  gr_mat_t recon;
  gr_mat_init(recon, gr_mat_nrows(ref, ctx), gr_mat_ncols(ref, ctx), ctx);

  status |= gr_struct_mat_reconstruct(recon, S, ctx);
  if (status != GR_SUCCESS || gr_mat_equal(recon, ref, ctx) == T_FALSE)
  {
    flint_printf("Reconstructed matrix does not match.\n");
    status = GR_TEST_FAIL;
  }

  slong current_rank = gr_mat_ncols(S->G, ctx);
  if (current_rank > max_expected_rank)
  {
    flint_printf("Rank after compression (%ld) exceeds expected max (%ld).\n", current_rank, max_expected_rank);
    status = GR_TEST_FAIL;
  }

  gr_mat_clear(recon, ctx);
  return status;
}

int test_compress_addition_toeplitz(int nb_iter)
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
    gr_mat_t A1, A2, Ref;

    gr_struct_mat_init(S1, n, m, T_TOEPLITZ, DISP_PLUS, ctx);
    gr_struct_mat_init(S2, n, m, T_TOEPLITZ, DISP_PLUS, ctx);
    gr_struct_mat_random(S1, state, ctx);
    gr_struct_mat_random(S2, state, ctx);

    gr_mat_init(A1, n, m, ctx);
    gr_mat_init(A2, n, m, ctx);
    gr_mat_init(Ref, n, m, ctx);

    gr_struct_mat_reconstruct(A1, S1, ctx);
    gr_struct_mat_reconstruct(A2, S2, ctx);
    status |= gr_mat_add(Ref, A1, A2, ctx); // Ref = S1 + S2

    status |= gr_struct_mat_add(S1, S2, ctx);
    if (status != GR_SUCCESS) break;

    status |= gr_struct_mat_compress(S1, ctx);
    if (status != GR_SUCCESS) break;

    status |= verify_compression(Ref, S1, 2, ctx);
    if (status != GR_SUCCESS) break;

    gr_struct_mat_clear(S1, ctx);
    gr_struct_mat_clear(S2, ctx);
    gr_mat_clear(A1, ctx);
    gr_mat_clear(A2, ctx);
    gr_mat_clear(Ref, ctx);

    gr_ctx_clear(ctx);
  }

  flint_rand_clear(state);
  return status;
}

int test_compress_multiplication_toeplitz(int nb_iter)
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
    slong p = n_randint(state, 100) + 10;

    gr_struct_mat_t S1, S2, S_prod;
    gr_mat_t A1, A2, Ref;

    gr_struct_mat_init(S1, n, m, T_TOEPLITZ, DISP_PLUS, ctx);
    gr_struct_mat_init(S2, m, p, T_TOEPLITZ, DISP_PLUS, ctx);
    gr_struct_mat_random(S1, state, ctx);
    gr_struct_mat_random(S2, state, ctx);

    gr_mat_init(A1, n, m, ctx);
    gr_mat_init(A2, m, p, ctx);
    gr_mat_init(Ref, n, p, ctx);

    gr_struct_mat_reconstruct(A1, S1, ctx);
    gr_struct_mat_reconstruct(A2, S2, ctx);
    status |= gr_mat_mul(Ref, A1, A2, ctx); // Ref = S1 * S2

    gr_struct_mat_init(S_prod, n, p, T_TOEPLITZ, DISP_PLUS, ctx);
    status |= gr_struct_mat_mul(S_prod, S1, S2, ctx);
    if (status != GR_SUCCESS) break;

    status |= gr_struct_mat_compress(S_prod, ctx);
    if (status != GR_SUCCESS) break;

    status |= verify_compression(Ref, S_prod, 4, ctx);
    if (status != GR_SUCCESS) break;

    gr_struct_mat_clear(S1, ctx);
    gr_struct_mat_clear(S2, ctx);
    gr_struct_mat_clear(S_prod, ctx);
    gr_mat_clear(A1, ctx);
    gr_mat_clear(A2, ctx);
    gr_mat_clear(Ref, ctx);

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

  if (strcmp("compress_addition_toeplitz", argv[1]) == 0) { ok = test_compress_addition_toeplitz(30); }
  else if (strcmp("compress_multiplication_toeplitz", argv[1]) == 0) { ok = test_compress_multiplication_toeplitz(30); }
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
