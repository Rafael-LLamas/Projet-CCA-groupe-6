#include "flint/flint.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "flint/gr_types.h"
#include "flint/nmod_mat.h"
#include "flint/nmod_poly.h"
#include "flint/ulong_extras.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "addition.h"
#include "displacement_matrices.h"
#include "matrix_aux.h"
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
  // Utilisation de CLOCK_MONOTONIC pour éviter les sauts de temps système
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
  if (!csv) return GR_UNABLE;

  // Tailles : Petite (128), Moyenne (1024), Grande (4096)
  slong sizes[] = {128, 1024, 4096};
  int num_sizes = 3;
  int iterations = 10; // Nombre de répétitions pour stats
  int error = GR_SUCCESS;
  fprintf(csv, "Size,Iteration,Gen_Time_ms,Flint_Dense_Time_ms\n");

  gr_ctx_t ctx;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

  printf("\n" ANSI_COLOR_BOLD ANSI_COLOR_MAGENTA
         "=== BENCHMARK MULTIPLICATION (Median/Avg over %d runs) ===" ANSI_COLOR_RESET "\n",
         iterations);
  printf("%-8s | %-12s | %-12s | %-12s\n", "Size", "Type", "Average (ms)", "Median (ms)");
  printf("---------|--------------|--------------|--------------\n");

  for (int s = 0; s < num_sizes; s++) {
    slong n = sizes[s];
    double gen_times[iterations];
    double flint_times[iterations];
    for (int i = 0; i < iterations; i++) {
      // 1. Initialiser matrices et générateurs aléatoires
      gr_mat_t C, A, B, G_a, G_b, H_a, H_b, G_c, H_c;

      gr_mat_init(A, n, n, ctx);
      gr_mat_init(B, n, n, ctx);
      gr_mat_init(C, n, n, ctx);
      error = random_toeplitz(A, n, n, state, ctx);
      error = random_toeplitz(B, n, n, state, ctx);
      error = gr_mat_G_H(G_a, H_a, A, ctx);
      error = gr_mat_G_H(G_b, H_b, B, ctx);
      // 2. Mesurer ma version
      double t1 = get_time_ms();
      error = gr_multiplication_generateur_deplacement_fast(G_c, H_c, G_a, H_a, G_b, H_b, ctx);
      gen_times[i] = get_time_ms() - t1;

      // 3. Mesurer FLINT
      double t3 = get_time_ms();
      error = gr_mat_mul(C, A, B, ctx);
      flint_times[i] = get_time_ms() - t3;

      // Log chaque itération dans le CSV
      fprintf(csv, "%ld,%d,%.4f,%.4f\n", n, i, gen_times[i], flint_times[i]);

      // Petit indicateur dynamique simple
      printf("\rTesting size %ld... [%d/%d]", n, i + 1, iterations);
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

    double avg_gen, med_gen, avg_flint, med_flint;
    compute_stats(gen_times, iterations, &avg_gen, &med_gen);
    compute_stats(flint_times, iterations, &avg_flint, &med_flint);

    // Affichage final pour cette taille
    printf("\r" ANSI_CLEAR_LINE); // Efface la ligne de chargement
    printf("%-8ld | " ANSI_COLOR_GREEN "Generators" ANSI_COLOR_RESET "   | %-12.3f | %-12.3f\n", n, avg_gen, med_gen);
    printf("%-8s | " ANSI_COLOR_YELLOW "Flint Dense" ANSI_COLOR_RESET "  | %-12.3f | %-12.3f\n", "", avg_flint,
           med_flint);
    printf("---------|--------------|--------------|--------------\n");

    // Écrire les stats finales pour cette taille dans le CSV
    fprintf(csv, "%ld,AVERAGE,%.4f,%.4f\n", n, avg_gen, avg_flint);
    fprintf(csv, "%ld,MEDIAN,%.4f,%.4f\n", n, med_gen, med_flint);
  }

  fclose(csv);
  gr_ctx_clear(ctx);
  flint_rand_clear(state);
  return error;
}
int benchmark_displacement() {
  FILE *csv = fopen("bench_displacement.csv", "w");
  if (!csv) return GR_UNABLE;

  slong sizes[] = {128, 512, 1024}; // On évite 4096 au début car la reconstruction est lente
  int num_sizes = 3;
  int iterations = 5;
  int error = GR_SUCCESS;

  fprintf(csv, "Size,Displacement_ms,Compression_GH_ms,Reconstruction_ms\n");

  gr_ctx_t ctx;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

  printf("\n" ANSI_COLOR_BOLD ANSI_COLOR_MAGENTA "=== BENCHMARK DISPLACEMENT ===" ANSI_COLOR_RESET "\n");
  printf("%-8s | %-12s | %-12s | %-12s\n", "Size", "Displace", "GH (Comp)", "Reconst");
  printf("---------|--------------|--------------|--------------\n");

  for (int s = 0; s < num_sizes; s++) {
    slong n = sizes[s];
    double t_disp = 0, t_gh = 0, t_rec = 0;

    for (int i = 0; i < iterations; i++) {
      gr_mat_t A, D, G, H, A_rec;
      gr_mat_init(A, n, n, ctx);
      gr_mat_init(D, n, n, ctx);
      gr_mat_init(G, 0, 0, ctx);
      gr_mat_init(H, 0, 0, ctx);
      gr_mat_init(A_rec, n, n, ctx);
      error = random_toeplitz(A, n, n, state, ctx);

      // 1. Temps Displacement
      double start = get_time_ms();
      error = gr_mat_displacement(D, A, ctx);
      t_disp += (get_time_ms() - start);

      // 2. Temps Compression (G et H)
      start = get_time_ms();
      error = gr_mat_G_H(G, H, A, ctx); // G et H sont init dans la fonction
      t_gh += (get_time_ms() - start);

      // 3. Temps Reconstruction
      start = get_time_ms();
      error = gr_mat_reconstruct_A_safe(A_rec, G, H, ctx);
      t_rec += (get_time_ms() - start);

      // Nettoyage
      gr_mat_clear(A, ctx);
      gr_mat_clear(D, ctx);
      gr_mat_clear(G, ctx);
      gr_mat_clear(H, ctx);
      gr_mat_clear(A_rec, ctx);
    }

    t_disp /= iterations;
    t_gh /= iterations;
    t_rec /= iterations;

    printf("%-8ld | %-12.3f | %-12.3f | %-12.3f\n", n, t_disp, t_gh, t_rec);
    fprintf(csv, "%ld,%.4f,%.4f,%.4f\n", n, t_disp, t_gh, t_rec);
  }

  fclose(csv);
  gr_ctx_clear(ctx);
  flint_rand_clear(state);
  return error;
}
int benchmark_addition() {
  FILE *csv = fopen("bench_addition.csv", "w");
  if (!csv) return GR_UNABLE;

  // Tailles : Petite (128), Moyenne (1024), Grande (4096)
  slong sizes[] = {128, 1024, 4096};
  int num_sizes = 3;
  int iterations = 10; // Nombre de répétitions pour stats
  int error = GR_SUCCESS;
  fprintf(csv, "Size,Iteration,Gen_Time_ms,Flint_Dense_Time_ms\n");

  gr_ctx_t ctx;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

  printf("\n" ANSI_COLOR_BOLD ANSI_COLOR_MAGENTA "=== BENCHMARK ADDITION (Median/Avg over %d runs) ===" ANSI_COLOR_RESET
         "\n",
         iterations);
  printf("%-8s | %-12s | %-12s | %-12s\n", "Size", "Type", "Average (ms)", "Median (ms)");
  printf("---------|--------------|--------------|--------------\n");

  for (int s = 0; s < num_sizes; s++) {
    slong n = sizes[s];
    double gen_times[iterations];
    double flint_times[iterations];
    for (int i = 0; i < iterations; i++) {
      // 1. Initialiser matrices et générateurs aléatoires
      gr_mat_t C, A, B, G_a, G_b, H_a, H_b, G_c, H_c;
      gr_mat_init(A, n, n, ctx);
      gr_mat_init(B, n, n, ctx);
      gr_mat_init(C, n, n, ctx);
      gr_mat_init(G_a, 0, 0, ctx);
      gr_mat_init(G_b, 0, 0, ctx);
      gr_mat_init(H_a, 0, 0, ctx);
      gr_mat_init(H_b, 0, 0, ctx);
      gr_mat_init(G_c, 0, 0, ctx);
      gr_mat_init(H_c, 0, 0, ctx);
      error = random_toeplitz(A, n, n, state, ctx);
      error = random_toeplitz(B, n, n, state, ctx);
      error = gr_mat_G_H(G_a, H_a, A, ctx);
      error = gr_mat_G_H(G_b, H_b, B, ctx);
      // 2. Mesurer ma version
      double t1 = get_time_ms();
      error = gr_mat_addition_generateur(G_c, H_c, G_a, H_a, G_b, H_b, ctx);
      gen_times[i] = get_time_ms() - t1;

      // 3. Mesurer FLINT
      double t3 = get_time_ms();
      error = gr_mat_add(C, A, B, ctx);
      flint_times[i] = get_time_ms() - t3;

      // Log chaque itération dans le CSV
      fprintf(csv, "%ld,%d,%.4f,%.4f\n", n, i, gen_times[i], flint_times[i]);

      // Petit indicateur dynamique simple
      printf("\rTesting size %ld... [%d/%d]", n, i + 1, iterations);
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

    double avg_gen, med_gen, avg_flint, med_flint;
    compute_stats(gen_times, iterations, &avg_gen, &med_gen);
    compute_stats(flint_times, iterations, &avg_flint, &med_flint);

    // Affichage final pour cette taille
    printf("\r" ANSI_CLEAR_LINE); // Efface la ligne de chargement
    printf("%-8ld | " ANSI_COLOR_GREEN "Generators" ANSI_COLOR_RESET "   | %-12.3f | %-12.3f\n", n, avg_gen, med_gen);
    printf("%-8s | " ANSI_COLOR_YELLOW "Flint Dense" ANSI_COLOR_RESET "  | %-12.3f | %-12.3f\n", "", avg_flint,
           med_flint);
    printf("---------|--------------|--------------|--------------\n");

    // Écrire les stats finales pour cette taille dans le CSV
    fprintf(csv, "%ld,AVERAGE,%.4f,%.4f\n", n, avg_gen, avg_flint);
    fprintf(csv, "%ld,MEDIAN,%.4f,%.4f\n", n, med_gen, med_flint);
  }

  fclose(csv);
  gr_ctx_clear(ctx);
  flint_rand_clear(state);
  return error;
};
void run_all_benchmarks() {
  printf("\033[H\033[J"); // Clear terminal
  printf(ANSI_COLOR_BOLD ANSI_COLOR_CYAN "=== GLOBAL BENCHMARK SUITE ===\n" ANSI_COLOR_RESET);
  benchmark_addition();
  benchmark_multiplication();
  benchmark_displacement();
  printf("\n" ANSI_COLOR_BOLD ANSI_COLOR_GREEN "Full report generated in .csv files.\n" ANSI_COLOR_RESET);
}
// Fonction d'usage
void usage(char *argv[]) {
  fprintf(stderr, ANSI_COLOR_BOLD ANSI_COLOR_CYAN "Usage: %s <command>\n" ANSI_COLOR_RESET, argv[0]);
  fprintf(stderr, "Available commands:\n");
  fprintf(stderr, "  " ANSI_COLOR_GREEN "benchmark" ANSI_COLOR_RESET " - Run all performance benchmarks\n");
  fprintf(stderr, "  " ANSI_COLOR_GREEN "benchmark_add" ANSI_COLOR_RESET "    - Run addition generator benchmark\n");
  fprintf(stderr,
          "  " ANSI_COLOR_GREEN "benchmark_mul" ANSI_COLOR_RESET "    - Run multiplication generator benchmark\n");
  fprintf(stderr,
          "  " ANSI_COLOR_GREEN "benchmark_displacement" ANSI_COLOR_RESET "    - Run displacement matrix benchmark\n");
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    usage(argv);
    return EXIT_FAILURE;
  }

  if (strcmp(argv[1], "benchmark") == 0) {
    printf(ANSI_COLOR_BOLD ANSI_COLOR_MAGENTA "Running all benchmarks...\n" ANSI_COLOR_RESET);
    run_all_benchmarks();
  } else if (strcmp(argv[1], "benchmark_add") == 0) {
    fprintf(stderr, ANSI_COLOR_BOLD ANSI_COLOR_BLUE "=> Running Addition Generator benchmark...\n" ANSI_COLOR_RESET);
    benchmark_addition();
  } else if (strcmp(argv[1], "benchmark_mul") == 0) {
    fprintf(stderr,
            ANSI_COLOR_BOLD ANSI_COLOR_BLUE "=> Running Multiplication Generator benchmark...\n" ANSI_COLOR_RESET);
    benchmark_multiplication();
  } else if (strcmp(argv[1], "benchmark_displacement") == 0) {
    fprintf(stderr, ANSI_COLOR_BOLD ANSI_COLOR_BLUE "=> Running Displacement Matrix benchmark...\n" ANSI_COLOR_RESET);
    benchmark_displacement();
  } else {
    fprintf(stderr, ANSI_COLOR_BOLD ANSI_COLOR_RED "Error:" ANSI_COLOR_RESET " Unknown command '%s'\n", argv[1]);
    usage(argv);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}