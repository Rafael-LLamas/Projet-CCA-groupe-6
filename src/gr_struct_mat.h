#ifndef GR_STRUCT_MAT_H
#define GR_STRUCT_MAT_H

#include "flint/gr.h"
#include "flint/gr_mat.h"
#include "flint/gr_types.h"

typedef enum
{
  T_TOEPLITZ,
  T_HANKEL,
  T_UNSURE
} structure_type_t;

typedef enum
{
  DISP_PLUS,
  DISP_MINUS
} disp_type_t;


typedef struct
{
  structure_type_t struct_t;
  disp_type_t disp_t;
  gr_mat_t G;
  gr_mat_t H;
} gr_struct_mat_struct;

typedef gr_struct_mat_struct gr_struct_mat_t[1];
typedef gr_struct_mat_struct *gr_struct_mat_ptr;
typedef const gr_struct_mat_struct *gr_struct_mat_srcptr;


// Basics
void gr_struct_mat_init(gr_struct_mat_t struct_mat, slong rows, slong cols, structure_type_t mat_struct, disp_type_t disp_t, gr_ctx_t ctx);
int gr_struct_mat_init_set(gr_struct_mat_t struct_mat, gr_mat_t mat, structure_type_t mat_struct, disp_type_t disp_t, gr_ctx_t ctx);
void gr_struct_mat_clear(gr_struct_mat_t mat, gr_ctx_t ctx);

slong gr_struct_mat_nrows(gr_struct_mat_t mat, gr_ctx_t ctx);
slong gr_struct_mat_ncols(gr_struct_mat_t mat, gr_ctx_t ctx);
slong gr_struct_mat_rank(gr_struct_mat_t mat, gr_ctx_t ctx);

void gr_struct_mat_print(gr_struct_mat_t mat, gr_ctx_t ctx);

int gr_struct_mat_reconstruct(gr_mat_t dense_mat, gr_struct_mat_t mat, gr_ctx_t ctx);

// Generator Compression
int gr_struct_mat_compress(gr_struct_mat_t mat, gr_ctx_t ctx);

// Addition
int gr_struct_mat_add(gr_struct_mat_t dest, gr_struct_mat_t mat, gr_ctx_t ctx);

#endif