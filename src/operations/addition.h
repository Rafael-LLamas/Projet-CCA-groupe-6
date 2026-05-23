#ifndef ADDITION_H
#define ADDITION_H

#include "flint/gr.h"

/**
 *
 * @brief the function will do an addition of A and B by concatenate their generator

 * @param[in] T Input the generator G of the matrix A
 * @param[in] U Input the generator H of the matrix A
 * @param[in] G Input the generator G of the matrix B
 * @param[in] H Input the generator H of the matrix B
 * @param[out] G_c Resulting G_a|G_b , the generator G of the matrix C
 * @param[out] H_c Resulting H_a|H_b , the generator H of the matrix C
 * @param[in] ctx The FLINT context object
 */
int gr_mat_addition_generateur(gr_mat_t T, gr_mat_t U, gr_mat_t G, gr_mat_t H, gr_mat_t G_c, gr_mat_t H_c,
                               gr_ctx_t ctx);

#endif
