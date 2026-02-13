#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "displacement_matrices.h"
#include "operations/multiplication.h"
#include "random_toeplitz.h"
#include "utility/matrix_aux.h"

typedef struct {
  const char *name;
  int (*func)(void);
} Command;

Command registry[] = {{"displacement", test_displacement_matrices},
                      {"random", test_random_toeplitz},
                      {"aux", test_matrix_aux},
                      {"multiplication", test_multiplication_generateurs},
                      {NULL, NULL}};

int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <test_name>\nAvailable tests:\n", argv[0]);
    for (int i = 0; registry[i].name != NULL; i++) { fprintf(stderr, "  - %s\n", registry[i].name); }
    return EXIT_FAILURE;
  }
  const char *test_target = argv[1];
  // run specific test
  for (int i = 0; registry[i].name != NULL; i++) {
    char filename[256];
    snprintf(filename, sizeof(filename), "output_%s.txt", test_target);
    if (freopen(filename, "w", stdout) == NULL) {
      perror("Failed to open test output file");
      return GR_TEST_FAIL;
    }
    int result = registry[i].func();
    if (result != 0) {
      fprintf(stderr, "FAILED: %s returned %d\n", test_target, result);
      return GR_TEST_FAIL;
    }

    fprintf(stderr, "PASSED: %s\n", test_target);
    return GR_SUCCESS;
  }
  fprintf(stderr, "Error: Test '%s' not found, maybe add it to registry[] dont forget cmake?\n", test_target);
  return EXIT_FAILURE;
}
