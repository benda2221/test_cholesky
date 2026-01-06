#include "test_helpers.h"
#include <string.h>

/* Matrix Operations */

void matrix_multiply_32x32(const double A[N][N], const double B[N][N], double C[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i][j] = 0.0;
            for (int k = 0; k < N; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void matrix_transpose_32x32(const double A[N][N], double At[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            At[i][j] = A[j][i];
        }
    }
}

bool matrix_compare_32x32(const double A[N][N], const double B[N][N], double tolerance) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (fabs(A[i][j] - B[i][j]) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

/* Matrix Generation */

void generate_identity_matrix(double A[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

void generate_random_spd_matrix(double A[N][N], unsigned int seed) {
    srand(seed);
    
    // Generate random lower triangular matrix L
    double L[N][N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (j > i) {
                L[i][j] = 0.0;  // Upper triangular part is zero
            } else {
                // Generate values in range [-1, 1] and add small positive bias
                L[i][j] = (2.0 * rand() / RAND_MAX - 1.0) * 0.5 + (i == j ? 0.5 : 0.0);
            }
        }
    }
    
    // Ensure diagonal is positive and not too small
    for (int i = 0; i < N; i++) {
        if (L[i][i] < 0.1) {
            L[i][i] = 0.1 + (double)rand() / RAND_MAX;
        }
    }
    
    // Compute A = L × L^T to get symmetric positive definite matrix
    double Lt[N][N];
    matrix_transpose_32x32(L, Lt);
    matrix_multiply_32x32(L, Lt, A);
}

/* Test Assertions */

bool assert_decomposition_correct(const double A[N][N], const double L[N][N],
                                  double tolerance, const char *test_name) {
    // Compute L × L^T
    double Lt[N][N];
    double reconstructed[N][N];
    
    matrix_transpose_32x32(L, Lt);
    matrix_multiply_32x32(L, Lt, reconstructed);
    
    // Compare with original A
    bool result = matrix_compare_32x32(A, reconstructed, tolerance);
    if (!result) {
        printf("FAIL: %s - Decomposition incorrect: L × L^T ≠ A\n", test_name);
        // Find maximum error
        double max_error = 0.0;
        int max_i = 0, max_j = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double error = fabs(A[i][j] - reconstructed[i][j]);
                if (error > max_error) {
                    max_error = error;
                    max_i = i;
                    max_j = j;
                }
            }
        }
        printf("  Maximum error: %e at [%d][%d]\n", max_error, max_i, max_j);
        printf("  A[%d][%d] = %e, (L×L^T)[%d][%d] = %e\n",
               max_i, max_j, A[max_i][max_j],
               max_i, max_j, reconstructed[max_i][max_j]);
    }
    return result;
}

bool assert_lower_triangular(const double L[N][N], double tolerance, const char *test_name) {
    bool result = true;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            if (fabs(L[i][j]) > tolerance) {
                if (result) {
                    printf("FAIL: %s - Matrix is not lower triangular\n", test_name);
                    result = false;
                }
                printf("  Upper triangular element [%d][%d] = %e (should be < %e)\n",
                       i, j, L[i][j], tolerance);
            }
        }
    }
    
    // Also check diagonal is positive
    for (int i = 0; i < N; i++) {
        if (L[i][i] <= 0.0) {
            if (result) {
                printf("FAIL: %s - Diagonal element is not positive\n", test_name);
                result = false;
            }
            printf("  Diagonal element [%d][%d] = %e (should be > 0)\n", i, i, L[i][i]);
        }
    }
    
    return result;
}

