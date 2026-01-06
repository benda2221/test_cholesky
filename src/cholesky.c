#include <math.h>
#include "cholesky.h"

#define N 32
#define EPSILON 1e-10

/**
 * @brief Performs Cholesky decomposition using Cholesky-Crout algorithm
 *
 * Algorithm:
 * For i = 0 to N-1:
 *   1. L[i][i] = sqrt(A[i][i] - sum(L[i][k]^2 for k < i))
 *   2. For j = i+1 to N-1:
 *      L[j][i] = (A[j][i] - sum(L[j][k] * L[i][k] for k < i)) / L[i][i]
 */
int cholesky_decompose_32x32(double A[32][32], double L[32][32]) {
    /* Step 1: Validate symmetry */
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            if (fabs(A[i][j] - A[j][i]) > EPSILON) {
                return 1; /* Matrix is not symmetric */
            }
        }
    }

    /* Step 2: Perform Cholesky decomposition */
    for (int i = 0; i < N; i++) {
        /* Compute diagonal element: L[i][i] = sqrt(A[i][i] - sum(L[i][k]^2)) */
        double sum = 0.0;
        for (int k = 0; k < i; k++) {
            sum += L[i][k] * L[i][k];
        }
        
        double diag_value = A[i][i] - sum;
        
        /* Check for positive definiteness */
        if (diag_value <= EPSILON) {
            return 2; /* Matrix is not positive definite */
        }
        
        L[i][i] = sqrt(diag_value);
        
        /* Compute off-diagonal elements in column i */
        for (int j = i + 1; j < N; j++) {
            sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += L[j][k] * L[i][k];
            }
            
            double numerator = A[j][i] - sum;
            
            /* Check for division by zero (should not happen if positive definite) */
            if (fabs(L[i][i]) < EPSILON) {
                return 3; /* Numerical error */
            }
            
            L[j][i] = numerator / L[i][i];
        }
        
        /* Zero upper triangular part for clarity and correctness
         * The algorithm only writes to lower triangular elements (L[j][i] where j >= i),
         * but we need upper triangular elements to be zero for L × L^T multiplication
         * to work correctly in verification and usage. */
        for (int j = 0; j < i; j++) {
            L[j][i] = 0.0;
        }
    }

    return 0; /* Success */
}

