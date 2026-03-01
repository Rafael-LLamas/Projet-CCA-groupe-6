#include "flint/flint.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "flint/ulong_extras.h"
#include <stdbool.h>
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

int size = -1;
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
  slong sizes[] = {128, 1024, 4096};
  int num_sizes = 3;
  if (!csv) return GR_UNABLE;

  if (size != -1) {
    sizes[0] = size / 2;
    sizes[1] = size;
    num_sizes = 2;
  }

  int iterations = 10;
  if (iteration != -1) { iterations = iteration; }

  fprintf(csv, "Size,Type,Iteration,My_Time_ms,Flint_Time_ms\n");

  gr_ctx_t ctx;
  flint_rand_t state;
  int error = 0;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

  printf("\n" ANSI_COLOR_BOLD ANSI_COLOR_MAGENTA "=== BENCHMARK MULTIPLICATION (Mat-Mat & Mat-Vec) ===" ANSI_COLOR_RESET
         "\n");

  for (int s = 0; s < num_sizes; s++) {
    slong n = sizes[s];
    double mat_mat_mine[iterations], mat_mat_flint[iterations];
    double mat_vec_mine[iterations], mat_vec_flint[iterations];

    for (int i = 0; i < iterations; i++) {
      gr_mat_t A, B, C, X, Res_V, G_a, H_a, G_b, H_b, G_c, H_c;

      // Init Matrices
      gr_mat_init(A, n, n, ctx);
      gr_mat_init(B, n, n, ctx);
      gr_mat_init(C, n, n, ctx);
      gr_mat_init(X, n, 1, ctx);
      gr_mat_init(Res_V, n, 1, ctx);

      // Init Générateurs
      gr_mat_init(G_a, 0, 0, ctx);
      gr_mat_init(H_a, 0, 0, ctx);
      gr_mat_init(G_b, 0, 0, ctx);
      gr_mat_init(H_b, 0, 0, ctx);
      gr_mat_init(G_c, 0, 0, ctx);
      gr_mat_init(H_c, 0, 0, ctx);

      // Données aléatoires
      error = random_toeplitz(A, n, n, state, ctx);
      error = random_toeplitz(B, n, n, state, ctx);
      error = gr_mat_G_H(G_a, H_a, A, ctx);
      error = gr_mat_G_H(G_b, H_b, B, ctx);
      // Remplir vecteur X
      for (slong r = 0; r < n; r++)
        error = gr_set(gr_mat_entry_ptr(X, r, 0, ctx), gr_mat_entry_srcptr(A, 0, r, ctx), ctx);

      // --- TEST 1 : MATRICE x MATRICE ---
      double t1 = get_time_ms();
      error = gr_multiplication_generateur_deplacement_fast(G_c, H_c, G_a, H_a, G_b, H_b, ctx);
      mat_mat_mine[i] = get_time_ms() - t1;

      double t2 = 0;
      if (flint) {
        t2 = get_time_ms();
        error = gr_mat_mul(C, A, B, ctx);
        mat_mat_flint[i] = get_time_ms() - t2;
      }
      mat_mat_flint[i] = t2 ? (get_time_ms() - t2) : 0;

      // --- TEST 2 : MATRICE x VECTEUR ---
      double t3 = get_time_ms();
      error = gr_mat_apply_struct_fast(Res_V, G_a, H_a, X, ctx);
      mat_vec_mine[i] = get_time_ms() - t3;

      double t4 = 0;
      if (flint) {
        t4 = get_time_ms();
        error = gr_mat_mul(Res_V, A, X, ctx); // Multiplication dense matrix-vector
        mat_vec_flint[i] = get_time_ms() - t4;
      }

      fprintf(csv, "%ld,MatMat,%d,%.4f,%.4f\n", n, i, mat_mat_mine[i], mat_mat_flint[i]);
      fprintf(csv, "%ld,MatVec,%d,%.4f,%.4f\n", n, i, mat_vec_mine[i], mat_vec_flint[i]);

      printf("\rSize %ld... [%d/%d]", n, i + 1, iterations);
      fflush(stdout);

      // Nettoyage complet
      gr_mat_clear(A, ctx);
      gr_mat_clear(B, ctx);
      gr_mat_clear(C, ctx);
      gr_mat_clear(X, ctx);
      gr_mat_clear(Res_V, ctx);
      gr_mat_clear(G_a, ctx);
      gr_mat_clear(H_a, ctx);
      gr_mat_clear(G_b, ctx);
      gr_mat_clear(H_b, ctx);
      gr_mat_clear(G_c, ctx);
      gr_mat_clear(H_c, ctx);
    }

    // Statistiques
    double avg_mm, med_mm, avg_mf = 0, med_mf = 0;
    double avg_vm, med_vm, avg_vf = 0, med_vf = 0;
    compute_stats(mat_mat_mine, iterations, &avg_mm, &med_mm);
    compute_stats(mat_vec_mine, iterations, &avg_vm, &med_vm);
    if (flint) {
      compute_stats(mat_mat_flint, iterations, &avg_mf, &med_mf);
      compute_stats(mat_vec_flint, iterations, &avg_vf, &med_vf);
    }

    // Affichage Table
    printf("\r" ANSI_CLEAR_LINE);
    printf(ANSI_COLOR_CYAN "Size %ld:\n" ANSI_COLOR_RESET, n);
    printf("  %-15s | %-12s | %-12s\n", "Operation", "Average (ms)", "Median (ms)");
    printf("  ----------------|--------------|--------------\n");
    printf("  " ANSI_COLOR_GREEN "Mat-Mat (Mine)" ANSI_COLOR_RESET "  | %-12.3f | %-12.3f\n", avg_mm, med_mm);
    if (flint)
      printf("  " ANSI_COLOR_YELLOW "Mat-Mat (Flint)" ANSI_COLOR_RESET " | %-12.3f | %-12.3f\n", avg_mf, med_mf);
    printf("  " ANSI_COLOR_BLUE "Mat-Vec (Mine)" ANSI_COLOR_RESET "  | %-12.3f | %-12.3f\n", avg_vm, med_vm);
    if (flint)
      printf("  " ANSI_COLOR_MAGENTA "Mat-Vec (Flint)" ANSI_COLOR_RESET " | %-12.3f | %-12.3f\n", avg_vf, med_vf);
    printf("  ----------------|--------------|--------------\n");
  }

  fclose(csv);
  gr_ctx_clear(ctx);
  flint_rand_clear(state);
  return error;
}
int benchmark_displacement() {
  FILE *csv = fopen("bench_displacement.csv", "w");
  slong sizes[] = {128, 512, 1024};
  int num_sizes = 3;
  if (!csv) return GR_UNABLE;
  if (size != -1) {
    sizes[0] = size / 2;
    sizes[1] = size;
    num_sizes = 2;
  }

  int iterations = 5;
  if (iteration != -1) { iterations = iteration; }
  int error = GR_SUCCESS;

  fprintf(csv, "Size,Displacement_ms,Compression_GH_ms,Reconstruction_ms\n");

  gr_ctx_t ctx;
  flint_rand_t state;
  flint_rand_init(state);
  flint_rand_set_seed(state, (ulong)time(NULL), (ulong)0x1234567890ABCDEF);
  gr_ctx_init_nmod(ctx, n_randprime(state, 64, 1));

  printf("\n" ANSI_COLOR_BOLD ANSI_COLOR_MAGENTA "=== BENCHMARK DISPLACEMENT ===" ANSI_COLOR_RESET "\n");
  printf("%-8s | %-12s | %-12s | %-12s\n", "Size", "Displace (ms)", "GH (ms)", "Reconst (ms)");
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
  slong sizes[] = {128, 1024, 4096};
  int num_sizes = 3;
  if (!csv) return GR_UNABLE;
  if (size != -1) {
    sizes[0] = size / 2;
    sizes[1] = size;
    num_sizes = 2;
  }

  int iterations = 10;
  if (iteration != -1) { iterations = iteration; }
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
      if (flint) {
        double t3 = get_time_ms();
        error = gr_mat_add(C, A, B, ctx);
        flint_times[i] = get_time_ms() - t3;
      }

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

    double avg_gen, med_gen;
    double avg_flint = 0.0, med_flint = 0.0;
    compute_stats(gen_times, iterations, &avg_gen, &med_gen);
    if (flint) { compute_stats(flint_times, iterations, &avg_flint, &med_flint); }

    // Affichage final pour cette taille
    printf("\r" ANSI_CLEAR_LINE); // Efface la ligne de chargement
    printf("%-8ld | " ANSI_COLOR_GREEN "Generators" ANSI_COLOR_RESET "   | %-12.3f | %-12.3f\n", n, avg_gen, med_gen);
    if (flint) {
      printf("%-8s | " ANSI_COLOR_YELLOW "Flint Dense" ANSI_COLOR_RESET "  | %-12.3f | %-12.3f\n", "", avg_flint,
             med_flint);
    }
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
  printf(ANSI_COLOR_BOLD ANSI_COLOR_CYAN "=== GLOBAL BENCHMARK ===\n" ANSI_COLOR_RESET);
  benchmark_addition();
  benchmark_multiplication();
  benchmark_displacement();
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
  fprintf(stderr,
          "  " ANSI_COLOR_GREEN "benchmark_displacement" ANSI_COLOR_RESET "    - Run displacement matrix benchmark\n");
  fprintf(stderr, "  " ANSI_COLOR_GREEN "-n (integer)" ANSI_COLOR_RESET
                  "    - Run matrix of size nxn with the integer (needs to be positive > 0)\n");
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
      size = atoi(argv[++i]);
      if (size <= 0) {
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
  } else if (strcmp(argv[1], "benchmark_displacement") == 0) {
    printf(ANSI_COLOR_BOLD ANSI_COLOR_BLUE "=> Running Displacement Matrix benchmark...\n" ANSI_COLOR_RESET);
    benchmark_displacement();
  } else {
    fprintf(stderr, ANSI_COLOR_BOLD ANSI_COLOR_RED "Error:" ANSI_COLOR_RESET " Unknown command '%s'\n", argv[1]);
    usage(argv);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}