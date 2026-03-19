#ifndef ADDITION_H
#define ADDITION_H

#include "flint/gr.h"

/**
 *
 * @brief the function will do an addition of A and B by concatenate their generator

 * @param[in] G_a Input the generator G of the matrix A
 * @param[in] H_a Input the generator H of the matrix A
 * @param[in] G_b Input the generator G of the matrix B
 * @param[in] H_b Input the generator H of the matrix B
 * @param[out] G_c Resulting G_a|G_b , the generator G of the matrix C
 * @param[out] H_c Resulting H_a|H_b , the generator H of the matrix C
 * @param[in] ctx The FLINT context object
 */
int gr_mat_addition_generateur(gr_mat_t G_a, gr_mat_t H_a, gr_mat_t G_b, gr_mat_t H_b, gr_mat_t G_c, gr_mat_t H_c,
                               gr_ctx_t ctx);

#endif
