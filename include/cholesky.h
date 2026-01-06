#ifndef CHOLESKY_H
#define CHOLESKY_H

/**
 * @file cholesky.h
 * @brief Cholesky decomposition for 32×32 positive definite matrices
 */

/**
 * @brief Performs Cholesky decomposition of a 32×32 symmetric positive definite matrix
 *
 * Decomposes matrix A into L × L^T where L is a lower triangular matrix.
 *
 * @param A Input: 32×32 symmetric positive definite matrix
 * @param L Output: 32×32 lower triangular matrix such that L × L^T = A
 * @return 0 on success, non-zero error code on failure:
 *         - 1: Matrix is not symmetric
 *         - 2: Matrix is not positive definite
 *         - 3: Numerical error during computation
 */
int cholesky_decompose_32x32(double A[32][32], double L[32][32]);

#endif /* CHOLESKY_H */

