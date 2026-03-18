#ifndef GR_MAT_RANDOM_TOEPLITZ_H
#define GR_MAT_RANDOM_TOEPLITZ_H

#include "flint/gr_mat.h"
/**
 * @brief the function will create a random quasi-toeplitz matrix of size n*m
 * @param[in] n Input row of A
 * @param[in] m Input collunm of A
 * @param[in] state The FLINT state for random value
 * @param[out] A Resulting a quasi-Toeplitz n*m
 * @param[in] ctx The FLINT context object
 */
int gr_mat_random_quasi_toepitz(gr_mat_t A, int n, int m, flint_rand_t state, gr_ctx_t ctx);
/**
 * @brief the function will create a random toeplitz matrix of size n*m
 * @param[in] n Input row of A
 * @param[in] m Input collunm of A
 * @param[in] state The FLINT state for random
 * @param[out] A Resulting a Toeplitz n*m
 * @param[in] ctx The FLINT context object
 */
int gr_mat_random_toeplitz(gr_mat_t A, int n, int m, flint_rand_t state, gr_ctx_t ctx);
/**
 * @brief the function will create a random quasi-toeplitz matrix of size n*m of rank nb_rand
 * @param[in] n Input row of A
 * @param[in] m Input collunm of A
 * @param[in] nb_rand the rank of A
 * @param[in] state The FLINT state for random
 * @param[out] A Resulting a quasi-Toeplitz n*m
 * @param[in] ctx The FLINT context object
 */
int gr_mat_Quasi_toeplitz_rank(gr_mat_t A, int n, int m, int nb_rand, flint_rand_t state, gr_ctx_t ctx);

#endif
