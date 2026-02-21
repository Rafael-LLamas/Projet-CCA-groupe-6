#include "displacement_matrices.h"
#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"
#include "matrix_aux.h"
#include "random_toeplitz.h"

/*
Arda's Notes:
It is mentioned in the book that we can use Strassen's Algorithm
to inverse a matrix.

Algorithm 10.1 Algorithme de type Strassen pour inverser une matrice
quasi-Toeplitz from [1].
*/

int gr_toeplitz_inverse(gr_mat_t G_D, gr_mat_t H_D, gr_mat_t G_A, gr_mat_t H_D, gr_ctx_t ctx) { return GR_SUCCESS; }
