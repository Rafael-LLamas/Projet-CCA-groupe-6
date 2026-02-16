#ifndef MATRIX_AUX_H
#define MATRIX_AUX_H

#include "flint/flint.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_poly.h"

#define GNMOD 1009

#define FLINT_CHECK(x)                                                                                                 \
  do {                                                                                                                 \
    int _status = (x);                                                                                                 \
    if (_status != GR_SUCCESS) return _status;                                                                         \
  } while (0)

int gr_mat_random_manual(gr_mat_t D, flint_rand_t state, ulong limit, gr_ctx_t ctx);

int gr_mat_lu_detach(gr_mat_t L, gr_mat_t U, gr_mat_t LU, gr_ctx_t ctx);

#endif