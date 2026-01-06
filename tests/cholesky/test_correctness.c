/**
 * @file test_correctness.c
 * @brief Correctness tests for Cholesky decomposition
 * 
 * Tests verify that the decomposition produces L such that L × L^T = A
 * and that L has the correct lower triangular structure.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../../include/cholesky.h"
#include "test_helpers.h"

/* Test counters */
static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

/* Test 1: Identity Matrix
 * Input: 32×32 identity matrix I
 * Expected: L = I (identity matrix)
 * Verification: L × L^T = I × I = I
 */
void test_identity_matrix(void) {
    tests_run++;
    printf("Test 1: Identity Matrix... ");
    
    double A[N][N];
    double L[N][N];
    
    generate_identity_matrix(A);
    
    int result = cholesky_decompose_32x32(A, L);
    
    if (result != 0) {
        printf("FAILED (error code %d)\n", result);
        tests_failed++;
        return;
    }
    
    // Verify L is identity
    if (!matrix_compare_32x32(A, L, TEST_TOLERANCE)) {
        printf("FAILED (L should equal identity)\n");
        tests_failed++;
        return;
    }
    
    // Verify decomposition correctness
    if (!assert_decomposition_correct(A, L, TEST_TOLERANCE, "Identity decomposition")) {
        tests_failed++;
        printf("FAILED\n");
        return;
    }
    
    // Verify lower triangular
    if (!assert_lower_triangular(L, TEST_TOLERANCE, "Identity lower triangular")) {
        tests_failed++;
        printf("FAILED\n");
        return;
    }
    
    tests_passed++;
    printf("PASSED\n");
}

/* Test 2: Diagonal Matrix
 * Input: Diagonal matrix D with positive diagonal elements
 * Expected: L is diagonal with sqrt of diagonal elements
 * Verification: L × L^T = D
 */
void test_diagonal_matrix(void) {
    tests_run++;
    printf("Test 2: Diagonal Matrix... ");
    
    double A[N][N];
    double L[N][N];
    
    // Create diagonal matrix with positive values
    memset(A, 0, sizeof(A));
    for (int i = 0; i < N; i++) {
        A[i][i] = (i + 1) * 2.0;  // Values: 2, 4, 6, ..., 64
    }
    
    int result = cholesky_decompose_32x32(A, L);
    
    if (result != 0) {
        printf("FAILED (error code %d)\n", result);
        tests_failed++;
        return;
    }
    
    // Verify L is diagonal with sqrt values
    bool diagonal_correct = true;
    for (int i = 0; i < N; i++) {
        double expected = sqrt(A[i][i]);
        if (fabs(L[i][i] - expected) > TEST_TOLERANCE) {
            printf("FAIL: Diagonal element [%d][%d] = %e, expected %e\n",
                   i, i, L[i][i], expected);
            diagonal_correct = false;
        }
        // Check off-diagonal is zero
        for (int j = 0; j < N; j++) {
            if (i != j && fabs(L[i][j]) > TEST_TOLERANCE) {
                printf("FAIL: Off-diagonal element [%d][%d] = %e (should be zero)\n",
                       i, j, L[i][j]);
                diagonal_correct = false;
            }
        }
    }
    
    if (!diagonal_correct) {
        tests_failed++;
        printf("FAILED\n");
        return;
    }
    
    // Verify decomposition correctness
    if (!assert_decomposition_correct(A, L, TEST_TOLERANCE, "Diagonal decomposition")) {
        tests_failed++;
        printf("FAILED\n");
        return;
    }
    
    tests_passed++;
    printf("PASSED\n");
}

/* Test 3: Random SPD Matrix
 * Input: Randomly generated symmetric positive definite matrix
 * Verification: L × L^T ≈ A (within numerical tolerance)
 * Verification: L is lower triangular
 */
void test_random_spd_matrix(void) {
    tests_run++;
    printf("Test 3: Random SPD Matrix... ");
    
    double A[N][N];
    double L[N][N];
    
    generate_random_spd_matrix(A, 12345);
    
    int result = cholesky_decompose_32x32(A, L);
    
    if (result != 0) {
        printf("FAILED (error code %d)\n", result);
        tests_failed++;
        return;
    }
    
    // Verify L is lower triangular
    if (!assert_lower_triangular(L, TEST_TOLERANCE, "Random SPD lower triangular")) {
        tests_failed++;
        printf("FAILED\n");
        return;
    }
    
    // Verify decomposition correctness
    if (!assert_decomposition_correct(A, L, TEST_TOLERANCE, "Random SPD decomposition")) {
        tests_failed++;
        printf("FAILED\n");
        return;
    }
    
    tests_passed++;
    printf("PASSED\n");
}

/* Test 4: Known Exact Decomposition
 * Input: Matrix with known exact Cholesky decomposition
 * Example: Simple 2×2 block repeated pattern
 * Verification: L × L^T = A exactly (within machine precision)
 */
void test_known_exact_decomposition(void) {
    tests_run++;
    printf("Test 4: Known Exact Decomposition... ");
    
    // Create a simple 2×2 block repeated pattern (scaled to 32×32)
    // For a 2×2 matrix [[a, b], [b, c]] with a>0 and ac-b^2>0,
    // the Cholesky decomposition is known
    
    double A[N][N];
    memset(A, 0, sizeof(A));
    
    // Create block diagonal with 2×2 blocks
    for (int block = 0; block < N/2; block++) {
        int i = block * 2;
        int j = block * 2 + 1;
        
        // Simple 2×2 SPD block: [[4, 1], [1, 4]]
        A[i][i] = 4.0;
        A[i][j] = 1.0;
        A[j][i] = 1.0;
        A[j][j] = 4.0;
    }
    
    double L[N][N];
    int result = cholesky_decompose_32x32(A, L);
    
    if (result != 0) {
        printf("FAILED (error code %d)\n", result);
        tests_failed++;
        return;
    }
    
    // Verify decomposition correctness with tight tolerance
    if (!assert_decomposition_correct(A, L, 1e-10, "Known exact decomposition")) {
        tests_failed++;
        printf("FAILED\n");
        return;
    }
    
    tests_passed++;
    printf("PASSED\n");
}

/* Test 5: Lower Triangular Property
 * Input: Random SPD matrix
 * Verification: Upper triangular part of L is zero
 * Verification: Diagonal elements are positive
 */
void test_lower_triangular_property(void) {
    tests_run++;
    printf("Test 5: Lower Triangular Property... ");
    
    double A[N][N];
    double L[N][N];
    
    generate_random_spd_matrix(A, 9999);
    
    int result = cholesky_decompose_32x32(A, L);
    
    if (result != 0) {
        printf("FAILED (error code %d)\n", result);
        tests_failed++;
        return;
    }
    
    // Verify lower triangular property
    if (!assert_lower_triangular(L, TEST_TOLERANCE, "Lower triangular property")) {
        tests_failed++;
        printf("FAILED\n");
        return;
    }
    
    tests_passed++;
    printf("PASSED\n");
}

/* Test 6: Multiple Random Matrices
 * Run decomposition on 10 different random SPD matrices
 * Verification: Correctness for each matrix
 */
void test_multiple_random_matrices(void) {
    tests_run++;
    printf("Test 6: Multiple Random Matrices (10 iterations)... ");
    
    int failed_count = 0;
    
    for (int iter = 0; iter < 10; iter++) {
        double A[N][N];
        double L[N][N];
        
        generate_random_spd_matrix(A, 1000 + iter);
        
        int result = cholesky_decompose_32x32(A, L);
        
        if (result != 0) {
            printf("FAIL: Iteration %d returned error code %d\n", iter, result);
            failed_count++;
            continue;
        }
        
        if (!assert_decomposition_correct(A, L, TEST_TOLERANCE, 
                                          "Multiple random matrices")) {
            failed_count++;
        }
    }
    
    if (failed_count > 0) {
        tests_failed++;
        printf("FAILED (%d/%d iterations failed)\n", failed_count, 10);
        return;
    }
    
    tests_passed++;
    printf("PASSED\n");
}

int main(void) {
    printf("=== Cholesky Decomposition Correctness Tests ===\n\n");
    
    test_identity_matrix();
    test_diagonal_matrix();
    test_random_spd_matrix();
    test_known_exact_decomposition();
    test_lower_triangular_property();
    test_multiple_random_matrices();
    
    printf("\n=== Test Summary ===\n");
    printf("Tests run: %d\n", tests_run);
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_failed);
    
    return (tests_failed == 0) ? 0 : 1;
}

