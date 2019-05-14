#include "dgemm.h"
#include <stdio.h>

const char * dgemm_desc = "DGEMM with loop indicies ordered for column major layout";

void square_dgemm(int N, double A[N][N], double B[N][N], double C[N][N])
{
	int i, j, k;
	for (k = 0; k < N; k++) {
		for (j = 0; j < N; j ++) {
			double tmp = B[k][j];
			for (i = 0; i < N; i++) {
				C[i][j] += A[i][k] * tmp;
			}
		}
	}
}
