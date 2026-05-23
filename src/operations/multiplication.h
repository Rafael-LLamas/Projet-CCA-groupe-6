#ifndef MULTIPLICATION_H
#define MULTIPLICATION_H

#include "flint/gr.h"

int test_multiplication_generateurs();

/**
 * @brief the function will do the multiplication of A*B by using their generators with W = Z*A*(Z^T)*G, V
 =(B^T)*U,a (resp. b) the last collunm of A (resp Z*(B^T))
 * @param[in] T Input the generator G of the matrix A
 * @param[in] U Input the generator H of the matrix A
 * @param[in] G Input the generator G of the matrix B
 * @param[in] H Input the generator H of the matrix B
 * @param[out] G_c Resulting T|W|alpha , the generator G of the matrix C
 * @param[out] H_c Resulting V|H|beta , the generator H of the matrix C
 * @param[in] ctx The FLINT context object
 */
int gr_mat_mul_generator(gr_mat_t G_c, gr_mat_t H_c, gr_mat_t T, gr_mat_t U, gr_mat_t G, gr_mat_t H,
                         gr_ctx_t ctx);
/**
 * @brief the function will multiply A*X with the generators of A and transforming X into a polynom
 * @param[in] G Input the generator G of the matrix A
 * @param[in] H Input the generator H of the matrix A
 * @param[in] V Input the vector V to multiply with A
 * @param[out] Res Resulting A*X
 * @param[in] ctx The FLINT context object
 */
int gr_mat_mul_vector(gr_mat_t Res, gr_mat_t G, gr_mat_t H, gr_mat_t V, gr_ctx_t ctx);
#endif
