#ifndef UTILITY_H
#define UTILITY_H

#include "gr_struct_mat.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"

int _gr_mat_lu_detach(gr_mat_t L, gr_mat_t U, gr_mat_t LU, gr_ctx_t ctx);

int _gr_mat_apply_Z(gr_mat_t Res, gr_mat_t M, gr_ctx_t ctx);
int _gr_mat_apply_Zt(gr_mat_t Res, gr_mat_t M, gr_ctx_t ctx);

int _gr_mat_displacement(gr_mat_t D, gr_mat_t A, disp_type_t type, gr_ctx_t ctx);

int _gr_mat_G_H(gr_mat_t G, gr_mat_t H, gr_mat_t A, disp_type_t type, gr_ctx_t ctx);

int _gr_struct_mat_guess_struct_t(structure_type_t *type, gr_mat_t mat, float acc_rate, gr_ctx_t ctx);

#endif