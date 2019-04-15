#include "dgemm.h"
#include <stdio.h>

const char * dgemm_desc = "DGEMM with loop indicies ordered for row major layout";

void square_dgemm(int N, double A[N][N], double B[N][N], double C[N][N])
{
	int i, j, k;
	for (k = 0; k < N; k++) {
		for (i = 0; i < N; i++) {
			double tmp = A[i][k];
			for (j = 0; j < N; j ++) {
				C[i][j] += tmp * B[k][j];
			}
		}
	}
}
