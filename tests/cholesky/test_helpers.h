#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 32
#define TEST_TOLERANCE 1e-9

/**
 * @file test_helpers.h
 * @brief Helper functions for Cholesky decomposition correctness tests
 */

/* Matrix Operations */

/**
 * @brief Multiply two 32×32 matrices: C = A × B
 * @param A First matrix (32×32)
 * @param B Second matrix (32×32)
 * @param C Output matrix (32×32) where C = A × B
 */
void matrix_multiply_32x32(const double A[N][N], const double B[N][N], double C[N][N]);

/**
 * @brief Compute transpose: At = A^T
 * @param A Input matrix (32×32)
 * @param At Output matrix (32×32) where At = A^T
 */
void matrix_transpose_32x32(const double A[N][N], double At[N][N]);

/**
 * @brief Compare two matrices with tolerance
 * @param A First matrix (32×32)
 * @param B Second matrix (32×32)
 * @param tolerance Maximum allowed difference
 * @return true if matrices are equal within tolerance, false otherwise
 */
bool matrix_compare_32x32(const double A[N][N], const double B[N][N], double tolerance);

/* Matrix Generation */

/**
 * @brief Generate 32×32 identity matrix
 * @param A Output matrix (32×32) initialized to identity
 */
void generate_identity_matrix(double A[N][N]);

/**
 * @brief Generate random symmetric positive definite matrix
 * @param A Output matrix (32×32) - will be SPD
 * @param seed Random seed for reproducibility
 * 
 * Algorithm: Generate random lower triangular matrix L, then compute A = L × L^T
 */
void generate_random_spd_matrix(double A[N][N], unsigned int seed);

/* Test Assertions */

/**
 * @brief Verify Cholesky decomposition correctness: L × L^T = A
 * @param A Original matrix (32×32)
 * @param L Lower triangular matrix from decomposition (32×32)
 * @param tolerance Maximum allowed difference
 * @param test_name Name of the test (for error messages)
 * @return true if decomposition is correct, false otherwise
 */
bool assert_decomposition_correct(const double A[N][N], const double L[N][N],
                                  double tolerance, const char *test_name);

/**
 * @brief Verify L is lower triangular (upper part is zero)
 * @param L Matrix to check (32×32)
 * @param tolerance Maximum allowed value in upper triangular part
 * @param test_name Name of the test (for error messages)
 * @return true if L is lower triangular, false otherwise
 */
bool assert_lower_triangular(const double L[N][N], double tolerance, const char *test_name);

#endif /* TEST_HELPERS_H */

