#ifndef MULTIPLICATION_H
#define MULTIPLICATION_H

#include "flint/gr.h"

int test_multiplication_generateurs();

/**
 * @brief the function will do the multiplication of A*B by using their generators with W = Z*A*(Z^T)*G_b, V
 =(B^T)*H_a,alpha (resp. beta) the last collunm of A (resp Z*(B^T))
 * @param[in] G_a Input the generator G of the matrix A
 * @param[in] H_a Input the generator H of the matrix A
 * @param[in] G_b Input the generator G of the matrix B
 * @param[in] H_b Input the generator H of the matrix B
 * @param[out] G_c Resulting G_a|W|alpha , the generator G of the matrix C
 * @param[out] H_c Resulting V|H_b|beta , the generator H of the matrix C
 * @param[in] ctx The FLINT context object
 */
int gr_mat_mul_generator(gr_mat_t G_c, gr_mat_t H_c, gr_mat_t G_a, gr_mat_t H_a, gr_mat_t G_b, gr_mat_t H_b,
                         gr_ctx_t ctx);
/**
 * @brief the function will multiply A*X with the generators of A and transforming X into a polynom
 * @param[in] G Input the generator G of the matrix A
 * @param[in] H Input the generator H of the matrix A
 * @param[in] X Input the vector X to multiply with A
 * @param[out] Res Resulting A*X
 * @param[in] ctx The FLINT context object
 */
int gr_mat_mul_vector(gr_mat_t Res, gr_mat_t G, gr_mat_t H, gr_mat_t X, gr_ctx_t ctx);
#endif
