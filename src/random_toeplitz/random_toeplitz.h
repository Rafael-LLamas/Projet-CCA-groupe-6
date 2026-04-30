#ifndef GR_MAT_RANDOM_TOEPLITZ_H
#define GR_MAT_RANDOM_TOEPLITZ_H

#include "flint/gr_mat.h"
/**
 * @brief the function will create a random quasi-toeplitz matrix of size n*m
 * @param[in] state The FLINT state for random value
 * @param[out] A Resulting a quasi-Toeplitz n*m
 * @param[in] ctx The FLINT context object
 */
int gr_mat_random_quasi_toepitz(gr_mat_t A, flint_rand_t state, gr_ctx_t ctx);
/**
 * @brief the function will create a random toeplitz matrix of size n*m
 * @param[in] state The FLINT state for random
 * @param[out] A Resulting a Toeplitz n*m
 * @param[in] ctx The FLINT context object
 */
int gr_mat_random_toeplitz(gr_mat_t A, flint_rand_t state, gr_ctx_t ctx);
/**
 * @brief the function will create a random quasi-toeplitz matrix of size n*m of rank nb_rand
 * @param[in] nb_rand the rank of A
 * @param[in] state The FLINT state for random
 * @param[out] A Resulting a quasi-Toeplitz n*m
 * @param[in] ctx The FLINT context object
 */
int gr_mat_quasi_toeplitz_rank(gr_mat_t A, int nb_rand, flint_rand_t state, gr_ctx_t ctx);
/**
 * @brief create directly random displacement generator for toeplitz matrix, G and H need to be init
 * @param[in] G one displacement generator of your future matrix
 * @param[in] H one displacement generator of your future matrix
 * @param[in] state The FLINT state for random
 * @param[in] ctx The FLINT context object
 */
int gr_mat_random_generator_toeplitz(gr_mat_t G, gr_mat_t H, flint_rand_t state, gr_ctx_t ctx);
#endif
