#include "dgemm.h"
#include <stdio.h>

const char * dgemm_desc = "DGEMM with unoptimized loop index order";

void square_dgemm(int N, double A[N][N], double B[N][N], double C[N][N])
{
	int i, j, k;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			double tmp = C[i][j];
			for (k = 0; k < N; k++) {
				tmp += A[i][k] * B[k][j];
			}
			C[i][j] = tmp;
		}
	}
}
