#include "dgemm.h"
#include <stdio.h>
#ifdef __INTEL_COMPILER
#include <mkl.h>
#else
#include <gsl/gsl_cblas.h>
#endif

const char * dgemm_desc = "DGEMM from system BLAS library";

void square_dgemm(int N, double A[N][N], double B[N][N], double C[N][N])
{
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0,
                    (double *)A, N, (double *)B, N, 1.0, (double *)C, N);
}
