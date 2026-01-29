#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "displacement_matrices.h"
#include "random_toeplitz.h"

struct Command {
  char *name;
  int (*func)(void);
};

// raphael add your functions here ---
struct Command registry[] = {
    {"displacement", test_displacement_matrices}, {"random", test_random_toeplitz}, {NULL, NULL}};

void usage() {
  fprintf(stderr, "Usage: ./main [test_name] ...\n");
  fprintf(stderr, "Available tests:\n");;
  for (int i = 0; registry[i].name != NULL; i++) { fprintf(stderr, "  - %s\n", registry[i].name); }
}

void execute_command(char *name, int (*func)(void)) {
  char filename[256];
  snprintf(filename, sizeof(filename), "output_%s.txt", name);
  if (freopen(filename, "w", stdout) == NULL) {
    perror("Failed to open output file");
    return;
  }
  if (func()) perror("A test failed.");
  fflush(stdout);
}

int run_selected(char *name) {
  for (int i = 0; registry[i].name != NULL; i++) {
    if (strcmp(name, registry[i].name) == 0) {
      execute_command(registry[i].name, registry[i].func);
      return 1;
    }
  }
  return 0;
}

int main(int argc, char *argv[]) {
  if (argc < 2) { // run all
    usage();
    for (int i = 0; registry[i].name != NULL; i++) execute_command(registry[i].name, registry[i].func);
    return EXIT_SUCCESS;
  }

  for (int i = 1; i < argc; i++) { // run selected
    if (!run_selected(argv[i])) fprintf(stderr, "Warning: Command '%s' not found.\n", argv[i]);
  }
  return EXIT_SUCCESS;
}