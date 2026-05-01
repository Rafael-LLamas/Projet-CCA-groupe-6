#include "flint/flint.h"
#include "flint/gr_mat.h"
#include "flint/ulong_extras.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "addition.h"
#include "compression.h"
#include "displacement_matrices.h"
#include "inverse_toeplitz.h"
#include "multiplication.h"
#include "random_toeplitz.h"

// --- GESTION DU TEMPS ---
#ifdef _WIN32
#include <windows.h>
double get_time_ms() {
  LARGE_INTEGER t, f;
  QueryPerformanceCounter(&t);
  QueryPerformanceFrequency(&f);
  return (double)t.QuadPart * 1000.0 / (double)f.QuadPart;
}
#else
#include <time.h>
double get_time_ms() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1000000.0;
}
#endif

// Codes ANSI pour un affichage joli
#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_BLUE "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_WHITE "\x1b[37m"
#define ANSI_COLOR_RESET "\x1b[0m"
#define ANSI_COLOR_BOLD "\x1b[1m"
#define ANSI_CURSOR_UP "\x1b[A"
#define ANSI_CLEAR_LINE "\x1b[2K"

int n = -1;
int m = -1;
int k = -1;
int rank = -1;
int iteration = -1;
bool flint = false;

void launch_external_terminal(int argc, char *argv[]) {
  char command[2048];
  char args_joined[512] = "";

  // On reconstruit les arguments SANS le "-out"
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-out") != 0) {
      strcat(args_joined, argv[i]);
      strcat(args_joined, " ");
    }
  }

#ifdef _WIN32
  // /Windows
  snprintf(command, sizeof(command), "start cmd /k \"%s %s\"", argv[0], args_joined);
#elif __APPLE__
  // MacOS
  snprintf(command, sizeof(command),
           "osascript -e 'tell application \"Terminal\" to do script \"cd \\\"$(pwd)\\\"; %s %s; echo; echo --- "
           "Termine ---; read -n 1 -s -p \\\"Appuyez sur une touche pour quitter...\\\"\"'",
           argv[0], args_joined);
#else
  // Linux
  snprintf(command, sizeof(command),
           "x-terminal-emulator -e bash -c \"./%s %s; echo; echo '--- Termine ---'; read -n 1 -s -r -p 'Appuyez sur "
           "une touche pour quitter...'\" 2>/dev/null &",
           argv[0], args_joined);
#endif

  printf(ANSI_COLOR_YELLOW "Lancement du benchmark dans le terminal externe...\n" ANSI_COLOR_RESET);
  int error = system(command);
  if (error != 0) { fprintf(stderr, "EROR ERROR ERRROR ERRRROR\n Try without -out"); }
}

int compare_doubles(const void *a, const void *b) {
  double arg1 = *(const double *)a;
  double arg2 = *(const double *)b;
  return (arg1 > arg2) - (arg1 < arg2);
}

void compute_stats(double *times, int count, double *avg, double *med) {
  double sum = 0;
  for (int i = 0; i < count; i++) sum += times[i];
  *avg = sum / count;

  qsort(times, count, sizeof(double), compare_doubles);

  if (count % 2 == 0)
    *med = (times[count / 2 - 1] + times[count / 2]) / 2.0;
  else
    *med = times[count / 2];
}

int benchmark_multiplication() {
  FILE *csv = fopen("bench_multiplication.csv", "w");
  slong sizesn[] = {128, 1024, 4096};
  slong sizesm[] = {128, 1024, 4096};
  slong sizesk[] = {128, 1024, 4096};
  int num_sizes = 3;
  int error = GR_SUCCESS;
  if (!csv) return GR_UNABLE;
  if (n != -1 || m != -1) {
    num_sizes = 2;
    slong base_n = (n != -1) ? n : m;
    slong base_m = (m != -1) ? m : n;
    sizesn[0] = base_n / 2;
    sizesn[1] = base_n;
    sizesm[0] = base_m / 2;
    sizesm[1] = base_m;
  }
  if (k != -1) {
    num_sizes = 2;
    sizesk[0] = k / 2;
    sizesk[1] = k;
  }

  int iterations = (iteration != -1) ? iteration : 10;
  fprintf(csv, "N,M,Type,Iteration,My_Time_ms,Flint_Time_ms\n");

  gr_ctx_t ctx;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

  printf("\n" ANSI_COLOR_BOLD ANSI_COLOR_MAGENTA "=== BENCHMARK MULTIPLICATION (Toeplitz & Random) ===" ANSI_COLOR_RESET
         "\n");

  for (int s = 0; s < num_sizes; s++) {
    slong ntemp = sizesn[s];
    slong mtemp = sizesm[s];
    slong ktemp = sizesk[s];

    double toe_toe_mine[iterations], toe_toe_flint[iterations];
    double toe_vec_mine[iterations], toe_vec_flint[iterations];

    for (int i = 0; i < iterations; i++) {
      gr_mat_t A, B, C, X, Res_V, Res_V2, G_a, H_a, G_b, H_b, G_c, H_c;

      gr_mat_init(A, ntemp, mtemp, ctx);
      gr_mat_init(B, mtemp, ktemp, ctx);
      gr_mat_init(C, ntemp, ktemp, ctx);
      gr_mat_init(X, mtemp, 1, ctx);
      gr_mat_init(Res_V, ntemp, 1, ctx);
      gr_mat_init(Res_V2, ntemp, 1, ctx);

      if (rank != -1) {
        error = gr_mat_quasi_toeplitz_rank(A, rank, state, ctx);
        error = gr_mat_quasi_toeplitz_rank(B, rank, state, ctx);
      } else {
        error = gr_mat_random_toeplitz(A, state, ctx);
        error = gr_mat_random_toeplitz(B, state, ctx);
      }

      error = gr_mat_G_H(G_a, H_a, A, DISP_PLUS, ctx);
      error = gr_mat_G_H(G_b, H_b, B, DISP_PLUS, ctx);
      for (slong r = 0; r < mtemp; r++)
        error = gr_set(gr_mat_entry_ptr(X, r, 0, ctx), gr_mat_entry_srcptr(A, 0, 0, ctx), ctx);

      double t1 = get_time_ms();
      error = gr_mat_mul_generator(G_c, H_c, G_a, H_a, G_b, H_b, ctx);
      toe_toe_mine[i] = get_time_ms() - t1;
      toe_toe_flint[i] = 0;
      if (flint) {
        double t2 = get_time_ms();
        error = gr_mat_mul(C, A, B, ctx);
        toe_toe_flint[i] = get_time_ms() - t2;
      }

      double t3 = get_time_ms();
      error = gr_mat_mul_vector(Res_V, G_a, H_a, X, ctx);
      toe_vec_mine[i] = get_time_ms() - t3;
      toe_vec_flint[i] = 0;
      if (flint) {
        double t4 = get_time_ms();
        error = gr_mat_mul(Res_V2, A, X, ctx);
        toe_vec_flint[i] = get_time_ms() - t4;
      }

      fprintf(csv, "%ld,%ld,ToeToe,%d,%.4f,%.4f\n", ntemp, mtemp, i, toe_toe_mine[i], toe_toe_flint[i]);
      fprintf(csv, "%ld,%ld,ToeVec,%d,%.4f,%.4f\n", ntemp, mtemp, i, toe_vec_mine[i], toe_vec_flint[i]);

      printf("\rSize %ldx%ld... [%d/%d]", ntemp, mtemp, i + 1, iterations);
      fflush(stdout);

      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_mat_clear(X, ctx);
      gr_mat_clear(Res_V, ctx);
      gr_mat_clear(Res_V2, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(G_b, ctx);
      gr_mat_clear(H_b, ctx);
      gr_mat_clear(G_c, ctx);
      gr_mat_clear(H_c, ctx);
    }

    double am_tt, mm_tt, af_tt = 0, mf_tt = 0, am_tv, mm_tv, af_tv = 0, mf_tv = 0;
    compute_stats(toe_toe_mine, iterations, &am_tt, &mm_tt);
    compute_stats(toe_vec_mine, iterations, &am_tv, &mm_tv);
    if (flint) {
      compute_stats(toe_toe_flint, iterations, &af_tt, &mf_tt);
      compute_stats(toe_vec_flint, iterations, &af_tv, &mf_tv);
    }

    printf("\r" ANSI_CLEAR_LINE);
    printf(ANSI_COLOR_CYAN "Size (%ldx%ld) X (%ldx%ld) :\n" ANSI_COLOR_RESET, ntemp, mtemp, mtemp, ktemp);
    printf("  %-15s | %-12s | %-12s\n", "Operation", "Average (ms)", "Median (ms)");
    printf("  ----------------|--------------|--------------\n");
    printf("  " ANSI_COLOR_GREEN "Toe-Toe (Mine)" ANSI_COLOR_RESET "  |  %-.3e   |  %-.3e \n", am_tt, mm_tt);
    if (flint)
      printf("  " ANSI_COLOR_YELLOW "Toe-Toe (Flint)" ANSI_COLOR_RESET " |  %-.3e   |  %-.3e \n", af_tt, mf_tt);
    printf("  " ANSI_COLOR_BLUE "Toe-Vec (Mine)" ANSI_COLOR_RESET "  |  %-.3e   |  %-.3e \n", am_tv, mm_tv);
    if (flint)
      printf("  " ANSI_COLOR_MAGENTA "Toe-Vec (Flint)" ANSI_COLOR_RESET " |  %-.3e   |  %-.3e \n", af_tv, mf_tv);
  }

  fclose(csv);
  gr_ctx_clear(ctx);
  flint_rand_clear(state);
  return error;
}
int benchmark_displacement() {
  FILE *csv = fopen("bench_auxiliary.csv", "w");
  slong sizes[] = {128, 512, 1024};
  int num_sizes = 3;
  if (!csv) return GR_UNABLE;
  if (n != -1 || m != -1) {
    slong base = (n > m) ? n : m;
    sizes[0] = base / 2;
    sizes[1] = base;
    num_sizes = 2;
  }

  int iterations = (iteration != -1) ? iteration : 5;
  int error = GR_SUCCESS;

  fprintf(csv, "Size,Displacement_ms,GH_ms,Compression_ms,Reconstruction_ms\n");

  gr_ctx_t ctx;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

  printf("\n" ANSI_COLOR_BOLD ANSI_COLOR_MAGENTA "=== BENCHMARK AUXILIARY ===" ANSI_COLOR_RESET "\n");

  for (int s = 0; s < num_sizes; s++) {
    slong cur_n = sizes[s];
    double t_disp = 0, t_gh = 0, t_rec = 0, t_comp = 0;
    ;
    double t_cur_disp = 0, t_cur_gh = 0, t_cur_rec = 0, t_cur_comp = 0;
    ;
    printf(ANSI_COLOR_CYAN "Size %ldx%ld :\n" ANSI_COLOR_RESET, cur_n, cur_n);
    printf(ANSI_COLOR_RED " %-12s " ANSI_COLOR_RESET "|" ANSI_COLOR_BLUE " %-12s " ANSI_COLOR_RESET "|" ANSI_COLOR_GREEN
                          " %-12s " ANSI_COLOR_RESET "|" ANSI_COLOR_YELLOW " %-12s\n" ANSI_COLOR_RESET,
           "Displace (ms)", "GH (ms)", "Compact (ms)", "Reconst (ms)");
    printf(" --------------|--------------|--------------|--------------\n");
    for (int i = 0; i < iterations; i++) {
      gr_mat_t A, D, G, H, A_rec, T, U, W, V;
      gr_mat_init(A, cur_n, cur_n, ctx);
      gr_mat_init(D, cur_n, cur_n, ctx);
      gr_mat_init(A_rec, cur_n, cur_n, ctx);

      if (rank != -1) {
        error = gr_mat_quasi_toeplitz_rank(A, rank, state, ctx);
      } else {
        error = gr_mat_random_toeplitz(A, state, ctx);
      }

      double start = get_time_ms();
      error = gr_mat_displacement(D, A, DISP_PLUS, ctx);
      t_cur_disp = (get_time_ms() - start);

      start = get_time_ms();
      error = gr_mat_G_H(G, H, A, DISP_PLUS, ctx);
      t_cur_gh = (get_time_ms() - start);
      error = gr_mat_init_set(T, G, ctx);
      error = gr_mat_init_set(U, H, ctx);
      error = gr_mat_zero(T, ctx);
      error = gr_mat_zero(U, ctx);
      error = gr_mat_addition_generateur(G, H, T, U, W, V, ctx);

      start = get_time_ms();
      error = gr_mat_generator_compress(W, V, ctx);
      t_cur_comp = (get_time_ms() - start);

      start = get_time_ms();
      error = gr_mat_reconstruct_A(A_rec, G, H, DISP_PLUS, ctx);
      t_cur_rec = (get_time_ms() - start);

      gr_mat_clear(A, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(G, ctx);
      gr_mat_clear(H, ctx);
      gr_mat_clear(T, ctx);
      gr_mat_clear(U, ctx);
      gr_mat_clear(W, ctx);
      gr_mat_clear(V, ctx);
      gr_mat_clear(A_rec, ctx);
      fprintf(csv, "%ld,%.4f,%.4f,%.4f,%.4f\n", cur_n, t_cur_disp, t_cur_gh, t_cur_comp, t_cur_rec);
      printf("\rSize %ldx%ld... [%d/%d]", cur_n, cur_n, i + 1, iterations);
      fflush(stdout);
      t_disp += t_cur_disp;
      t_gh += t_cur_gh;
      t_rec += t_cur_rec;
      t_comp += t_cur_comp;
    }

    t_disp /= iterations;
    t_gh /= iterations;
    t_rec /= iterations;
    t_comp /= iterations;
    printf("\r" ANSI_CLEAR_LINE);
    printf("  %-.3e    |  %-.3e   |  %-.3e   |  %-.3e \n", t_disp, t_gh, t_comp, t_rec);
  }

  fclose(csv);
  gr_ctx_clear(ctx);
  flint_rand_clear(state);
  return error;
}

int benchmark_addition() {
  FILE *csv = fopen("bench_addition.csv", "w");
  slong sizesn[] = {128, 1024, 4096};
  slong sizesm[] = {128, 1024, 4096};
  int num_sizes = 3;
  if (!csv) return GR_UNABLE;

  if (n != -1 || m != -1) {
    num_sizes = 2;
    slong base_n = (n != -1) ? n : m;
    slong base_m = (m != -1) ? m : n;
    sizesn[0] = base_n / 2;
    sizesn[1] = base_n;
    sizesm[0] = base_m / 2;
    sizesm[1] = base_m;
  }

  int iterations = (iteration != -1) ? iteration : 10;
  int error = GR_SUCCESS;
  fprintf(csv, "N,M,Iteration,Gen_Time_ms,Flint_Dense_Time_ms\n");

  gr_ctx_t ctx;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

  printf("\n" ANSI_COLOR_BOLD ANSI_COLOR_MAGENTA "=== BENCHMARK ADDITION (%d runs) ===" ANSI_COLOR_RESET "\n",
         iterations);

  for (int s = 0; s < num_sizes; s++) {
    slong cur_n = sizesn[s];
    slong cur_m = sizesm[s];
    double gen_times[iterations], flint_times[iterations];
    printf(ANSI_COLOR_CYAN "Size %ldx%ld :\n" ANSI_COLOR_RESET, cur_n, cur_m);
    printf(" %-12s | %-12s | %-12s\n", "Operation", "Average (ms)", "Median (ms)");
    printf("--------------|--------------|--------------\n");

    for (int i = 0; i < iterations; i++) {
      gr_mat_t C, A, B, G_a, G_b, H_a, H_b, G_c, H_c;
      gr_mat_init(A, cur_n, cur_m, ctx);
      gr_mat_init(B, cur_n, cur_m, ctx);
      gr_mat_init(C, cur_n, cur_m, ctx);
      if (rank != -1) {
        error = gr_mat_quasi_toeplitz_rank(A, rank, state, ctx);
        error = gr_mat_quasi_toeplitz_rank(B, rank, state, ctx);

      } else {
        error = gr_mat_random_toeplitz(A, state, ctx);
        error = gr_mat_random_toeplitz(B, state, ctx);
      }
      error = gr_mat_G_H(G_a, H_a, A, DISP_PLUS, ctx);
      error = gr_mat_G_H(G_b, H_b, B, DISP_PLUS, ctx);
      double t1 = get_time_ms();
      error = gr_mat_addition_generateur(G_a, H_a, G_b, H_b, G_c, H_c, ctx);
      gen_times[i] = get_time_ms() - t1;

      flint_times[i] = 0;
      if (flint) {
        double t3 = get_time_ms();
        error = gr_mat_add(C, A, B, ctx);
        flint_times[i] = get_time_ms() - t3;
      }

      fprintf(csv, "%ld,%ld,%d,%.4f,%.4f\n", cur_n, cur_m, i, gen_times[i], flint_times[i]);
      printf("\rSize %ldx%ld... [%d/%d]", cur_n, cur_m, i + 1, iterations);
      fflush(stdout);

      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(G_b, ctx);
      gr_mat_clear(H_b, ctx);
      gr_mat_clear(G_c, ctx);
      gr_mat_clear(H_c, ctx);
    }

    double avg_gen, med_gen, avg_flint = 0, med_flint = 0;
    compute_stats(gen_times, iterations, &avg_gen, &med_gen);
    if (flint) compute_stats(flint_times, iterations, &avg_flint, &med_flint);

    printf("\r" ANSI_CLEAR_LINE);
    printf(ANSI_COLOR_GREEN " Generators" ANSI_COLOR_RESET "   |  %-.3e   |  %-.3e \n", avg_gen, med_gen);
    if (flint)
      printf(ANSI_COLOR_YELLOW " Flint Dense" ANSI_COLOR_RESET "  |  %-.3e   |  %-.3e \n", avg_flint, med_flint);
  }

  fclose(csv);
  gr_ctx_clear(ctx);
  flint_rand_clear(state);
  return error;
}

int benchmark_inversion() {
  FILE *csv = fopen("bench_inversion.csv", "w");
  slong sizesn[] = {512, 1024, 2048};
  int num_sizes = 3;
  if (!csv) return GR_UNABLE;

  if (n != -1) {
    num_sizes = 2;
    sizesn[0] = n / 2;
    sizesn[1] = n;
  }

  int iterations = (iteration != -1) ? iteration : 10;
  int error = GR_SUCCESS;
  fprintf(csv, "N,M,Iteration,Gen_Time_ms,Flint_Time_ms\n");

  gr_ctx_t ctx;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

  printf("\n" ANSI_COLOR_BOLD ANSI_COLOR_MAGENTA "=== BENCHMARK INVERSSION (%d runs) ===" ANSI_COLOR_RESET "\n",
         iterations);

  for (int s = 0; s < num_sizes; s++) {
    slong cur_n = sizesn[s];
    double gen_times[iterations], flint_times[iterations];
    printf(ANSI_COLOR_CYAN "Size %ldx%ld :\n" ANSI_COLOR_RESET, cur_n, cur_n);
    printf(" %-12s| %-12s | %-12s\n", "Operation", "Average (ms)", "Median (ms)");
    printf("-------------|--------------|--------------\n");

    for (int i = 0; i < iterations; i++) {
      gr_mat_t A, B, G_a, G_b, H_a, H_b;
      gr_ptr det;
      gr_mat_init(A, cur_n, cur_n, ctx);
      gr_mat_init(B, cur_n, cur_n, ctx);
      GR_TMP_INIT(det, ctx);

      if (rank != -1) {
        error = gr_mat_quasi_toeplitz_rank(A, rank, state, ctx);
      } else {
        error = gr_mat_random_toeplitz(A, state, ctx);
      }

      error = gr_mat_G_H(G_a, H_a, A, DISP_PLUS, ctx);
      double t1 = get_time_ms();
      error = gr_mat_inverse_toeplitz(G_b, H_b, G_a, H_a, ctx);
      if (error) {
        GR_TMP_CLEAR(det, ctx);
        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(G_a, ctx);
        gr_mat_clear(H_a, ctx);
        gr_mat_clear(G_b, ctx);
        gr_mat_clear(H_b, ctx);
        continue;
      }
      gen_times[i] = get_time_ms() - t1;

      flint_times[i] = 0;
      if (flint) {
        double t3 = get_time_ms();
        error = gr_mat_inv(B, A, ctx);
        flint_times[i] = get_time_ms() - t3;
      }

      fprintf(csv, "%ld,%ld,%d,%.4f,%.4f\n", cur_n, cur_n, i, gen_times[i], flint_times[i]);
      printf("\rSize %ldx%ld... [%d/%d]", cur_n, cur_n, i + 1, iterations);
      fflush(stdout);
      GR_TMP_CLEAR(det, ctx);
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(G_b, ctx);
      gr_mat_clear(H_b, ctx);
    }

    double avg_gen, med_gen, avg_flint = 0, med_flint = 0;
    compute_stats(gen_times, iterations, &avg_gen, &med_gen);
    if (flint) compute_stats(flint_times, iterations, &avg_flint, &med_flint);

    printf("\r" ANSI_CLEAR_LINE);
    printf(ANSI_COLOR_GREEN " Generators " ANSI_COLOR_RESET " |  %-.3e   |  %-.3e \n", avg_gen, med_gen);
    if (flint)
      printf(ANSI_COLOR_YELLOW " Flint Dense" ANSI_COLOR_RESET " |  %-.3e   |  %-.3e \n", "", avg_flint, med_flint);
  }

  fclose(csv);
  gr_ctx_clear(ctx);
  flint_rand_clear(state);
  return error;
}
int benchmark_random() {
  FILE *csv = fopen("bench_random.csv", "w");
  slong sizesn[] = {128, 512, 1024};
  slong sizesm[] = {128, 512, 1024};
  int num_sizes = 3;
  if (!csv) return GR_UNABLE;

  if (n != -1 || m != -1) {
    num_sizes = 2;
    slong base_n = (n != -1) ? n : m;
    slong base_m = (m != -1) ? m : n;
    sizesn[0] = base_n / 2;
    sizesn[1] = base_n;
    sizesm[0] = base_m / 2;
    sizesm[1] = base_m;
  }

  int iterations = (iteration != -1) ? iteration : 10;
  if (rank == -1) { rank = 3; }
  int error = GR_SUCCESS;
  fprintf(csv, "N,M,Iteration,Gen_Toeplitz_Time_ms,Matrice_Toeplitz_Dense_time_ms,Matrice_Quasi_Toeplitz_Dense_time_ms,"
               "Matrice_Quasi_Toeplitz_rank_Dense_time_ms\n");

  gr_ctx_t ctx;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

  printf("\n" ANSI_COLOR_BOLD ANSI_COLOR_MAGENTA "=== BENCHMARK RANDOM (%d runs) ===" ANSI_COLOR_RESET "\n",
         iterations);

  for (int s = 0; s < num_sizes; s++) {
    slong cur_n = sizesn[s];
    slong cur_m = sizesm[s];
    double gen_Toe[iterations], Toe_dense[iterations], quasi_Toe_dense[iterations], quasi_Toe_rank_dense[iterations];
    printf(ANSI_COLOR_CYAN "Size %ldx%ld :\n" ANSI_COLOR_RESET, cur_n, cur_m);
    printf(" %-12s   | %-12s | %-12s\n", "Operation", "Average (ms)", "Median (ms)");
    printf("----------------|--------------|--------------\n");

    for (int i = 0; i < iterations; i++) {
      gr_mat_t C, A, B, G, H, T, U;
      gr_mat_init(A, cur_n, cur_m, ctx);
      gr_mat_init(B, cur_n, cur_m, ctx);
      gr_mat_init(C, cur_n, cur_m, ctx);
      gr_mat_init(G, cur_n, 2, ctx);
      gr_mat_init(H, cur_m, 2, ctx);

      double t1 = get_time_ms();
      error = gr_mat_random_toeplitz(A, state, ctx);
      Toe_dense[i] = get_time_ms() - t1;
      double t2 = get_time_ms();
      error = gr_mat_random_quasi_toepitz(B, state, ctx);
      quasi_Toe_dense[i] = get_time_ms() - t2;
      double t3 = get_time_ms();
      error = gr_mat_quasi_toeplitz_rank(C, rank, state, ctx);
      quasi_Toe_rank_dense[i] = get_time_ms() - t3;
      double t4 = get_time_ms();
      error = gr_mat_random_generator_toeplitz(G, H, state, ctx);
      gen_Toe[i] = get_time_ms() - t4;

      fprintf(csv, "%ld,%ld,%d,%.4f,%.4f,%.4f,%.4f\n", cur_n, cur_m, i, gen_Toe[i], Toe_dense[i], quasi_Toe_dense[i],
              quasi_Toe_rank_dense[i]);
      printf("\rSize %ldx%ld... [%d/%d]", cur_n, cur_m, i + 1, iterations);
      fflush(stdout);

      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_mat_clear(G, ctx);
      gr_mat_clear(H, ctx);
    }

    double avg_gen_toe, med_gen_toe, avg_toe_dense, med_toe_dense, avg_quasi_dense, med_quasi_dense, avg_rank, med_rank;
    compute_stats(gen_Toe, iterations, &avg_gen_toe, &med_gen_toe);
    compute_stats(Toe_dense, iterations, &avg_toe_dense, &med_toe_dense);
    compute_stats(quasi_Toe_dense, iterations, &avg_quasi_dense, &med_quasi_dense);
    compute_stats(quasi_Toe_rank_dense, iterations, &avg_rank, &med_rank);

    printf("\r" ANSI_CLEAR_LINE);
    printf(ANSI_COLOR_GREEN " Generators Toe" ANSI_COLOR_RESET " |  %-.3e   |  %-.3e \n", avg_gen_toe, med_gen_toe);
    printf(ANSI_COLOR_BLUE " Matrix Toe" ANSI_COLOR_RESET "     |  %-.3e   |  %-.3e \n", avg_toe_dense, med_toe_dense);
    printf(ANSI_COLOR_MAGENTA " Matrix Quasi" ANSI_COLOR_RESET "   |  %-.3e   |  %-.3e \n", avg_quasi_dense,
           med_quasi_dense);
    printf(ANSI_COLOR_YELLOW " Matrix Rank %d" ANSI_COLOR_RESET "  |  %-.3e   |  %-.3e \n", rank, avg_rank, med_rank);
  }

  fclose(csv);
  gr_ctx_clear(ctx);
  flint_rand_clear(state);
  return error;
}
void run_all_benchmarks() {
  printf("\033[H\033[J"); // Clear terminal
  printf(ANSI_COLOR_BOLD ANSI_COLOR_CYAN "=== GLOBAL BENCHMARK ===\n" ANSI_COLOR_RESET);
  benchmark_addition();
  benchmark_multiplication();
  benchmark_displacement();
  benchmark_inversion();
  benchmark_random();
  printf("\n" ANSI_COLOR_BOLD ANSI_COLOR_GREEN "Full report generated in .csv files.\n" ANSI_COLOR_RESET);
}
// Fonction d'usage
void usage(char *argv[]) {
  fprintf(stderr, ANSI_COLOR_BOLD ANSI_COLOR_CYAN "Usage: %s <command> <options>\n" ANSI_COLOR_RESET, argv[0]);
  fprintf(stderr, "Available commands:\n");
  fprintf(stderr, "  " ANSI_COLOR_GREEN "benchmark" ANSI_COLOR_RESET " - Run all performance benchmarks\n");
  fprintf(stderr, "  " ANSI_COLOR_GREEN "benchmark_add" ANSI_COLOR_RESET "    - Run addition generator benchmark\n");
  fprintf(stderr,
          "  " ANSI_COLOR_GREEN "benchmark_mul" ANSI_COLOR_RESET "    - Run multiplication generator benchmark\n");
  fprintf(stderr, "  " ANSI_COLOR_GREEN "benchmark_inv" ANSI_COLOR_RESET "    - Run inversion generator benchmark\n");
  fprintf(stderr, "  " ANSI_COLOR_GREEN "benchmark_aux" ANSI_COLOR_RESET "    - Run auxiliary matrix benchmark\n");
  fprintf(stderr, "  " ANSI_COLOR_GREEN "benchmark_rand" ANSI_COLOR_RESET "    - Run random matrix benchmark\n");
  fprintf(stderr, "  " ANSI_COLOR_GREEN "-n (integer)" ANSI_COLOR_RESET
                  "    - Run matrix of nxm or nxn (if m not define) with the integer (needs to be positive > 0)\n");
  fprintf(stderr, "  " ANSI_COLOR_GREEN "-k (integer)" ANSI_COLOR_RESET
                  "    - Run matrix of mxk (if m not define) with the integer (needs to be positive > 0)\n");
  fprintf(stderr, "  " ANSI_COLOR_GREEN "-m (integer)" ANSI_COLOR_RESET
                  "    - Run matrix of nxm or mxm (if n not define) with the integer (needs to be positive > 0)\n");
  fprintf(stderr, "  " ANSI_COLOR_GREEN "-i (integer)" ANSI_COLOR_RESET
                  "    - Run n iterations with the integer (needs to be positive > 0)\n");
  fprintf(stderr, "  " ANSI_COLOR_GREEN "-r (integer)" ANSI_COLOR_RESET
                  "    - Run matrix of rank n with the integer (needs to be positive > 0)\n");
  fprintf(stderr,
          "  " ANSI_COLOR_GREEN "-out" ANSI_COLOR_RESET "    - Run a new terminal (only work with Debian like's )\n");
  fprintf(stderr, "  " ANSI_COLOR_GREEN "-flint" ANSI_COLOR_RESET
                  "    - Compare multiplication and addition with flint times\n");
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    usage(argv);
    return EXIT_FAILURE;
  }

  bool external = false;

  // --- PARSING DES ARGUMENTS ---
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-out") == 0) {
      external = true;
    } else if (strcmp(argv[i], "-flint") == 0) {
      flint = true;
    } else if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
      n = atoi(argv[++i]);
      if (n <= 0) {
        fprintf(stderr, "Error: -n needs to be > 0\n");
        return EXIT_FAILURE;
      }
    } else if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
      iteration = atoi(argv[++i]);
      if (iteration <= 0) {
        fprintf(stderr, "Error: -i needs to be > 0\n");
        return EXIT_FAILURE;
      }
    } else if (strcmp(argv[i], "-r") == 0 && i + 1 < argc) {
      rank = atoi(argv[++i]);
      if (rank <= 0) {
        fprintf(stderr, "Error: -r needs to be > 0\n");
        return EXIT_FAILURE;
      }
    } else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
      m = atoi(argv[++i]);
      if (m <= 0) {
        fprintf(stderr, "Error: -m needs to be > 0\n");
        return EXIT_FAILURE;
      }
    } else if (strcmp(argv[i], "-k") == 0 && i + 1 < argc) {
      k = atoi(argv[++i]);
      if (k <= 0) {
        fprintf(stderr, "Error: -k needs to be > 0\n");
        return EXIT_FAILURE;
      }
    }
  }

  // Gestion du terminal externe
  if (external) {
    launch_external_terminal(argc, argv);
    printf(ANSI_COLOR_GREEN "Le terminal externe a été lancé.\n" ANSI_COLOR_RESET);
    printf("Appuyez sur Entrée pour fermer ce terminal...\n");
    getchar();
    return EXIT_SUCCESS;
  }

  // --- EXECUTION DES COMMANDES ---
  // La commande est toujours le premier argument non-optionnel (argv[1])
  if (strcmp(argv[1], "benchmark") == 0) {
    run_all_benchmarks();
  } else if (strcmp(argv[1], "benchmark_add") == 0) {
    printf(ANSI_COLOR_BOLD ANSI_COLOR_BLUE "=> Running Addition Generator benchmark...\n" ANSI_COLOR_RESET);
    benchmark_addition();
  } else if (strcmp(argv[1], "benchmark_mul") == 0) {
    printf(ANSI_COLOR_BOLD ANSI_COLOR_BLUE "=> Running Multiplication Generator benchmark...\n" ANSI_COLOR_RESET);
    benchmark_multiplication();
  } else if (strcmp(argv[1], "benchmark_inv") == 0) {
    printf(ANSI_COLOR_BOLD ANSI_COLOR_BLUE "=> Running Inversion Generator benchmark...\n" ANSI_COLOR_RESET);
    benchmark_inversion();
  } else if (strcmp(argv[1], "benchmark_aux") == 0) {
    printf(ANSI_COLOR_BOLD ANSI_COLOR_BLUE "=> Running Auxiliary Matrix benchmark...\n" ANSI_COLOR_RESET);
    benchmark_displacement();
  } else if (strcmp(argv[1], "benchmark_rand") == 0) {
    printf(ANSI_COLOR_BOLD ANSI_COLOR_BLUE "=> Running Random Matrix benchmark...\n" ANSI_COLOR_RESET);
    benchmark_random();
  } else {
    fprintf(stderr, ANSI_COLOR_BOLD ANSI_COLOR_RED "Error:" ANSI_COLOR_RESET " Unknown command '%s'\n", argv[1]);
    usage(argv);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
