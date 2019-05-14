#include "dgemm.h"
#include <stdio.h>

/* Cache blocking parameters */
#ifndef BI
#define BI 96
#endif

#ifndef BJ
#define BJ 96
#endif 

#ifndef BK
#define BK 32
#endif

#define xstr(s) str(s)
#define str(s) #s

#define CBLOCK_DESCR BI * BJ * BK

const char * dgemm_desc = "DGEMM cache blocked " xstr(CBLOCK_DESCR);

inline int min(int i, int j)
{
	return i < j ? i : j;
}

/*****
 * General, mildly optimized versions of dgebb and dgepb
 * that can be used in the general case. An optimized implementation
 * would register block and vectorize dgebb for specific block sizes.
 */

/* Multiply a BI*BK block of A times a BK*BJ block of B 
 * into an BI * BJ block of C. In the common case, BI = BJ = BK */
void dgebb(int bi, int bj, int bk,
	   int Astride, double A[][Astride],  
	   int Bstride, double B[][Bstride], 
	   int Cstride, double C[][Cstride])
{
	int i, j, k;
	for (i = 0; i < bi; i++) {
		for (k = 0; k < bk; k++) {
			double tmp = A[i][k];
			for (j = 0; j < bj; j++) {
				C[i][j] += tmp * B[k][j];
			}
		}
	}
}

/* Multiply a NxBK panel of A times a BKxBJ block of B
 * into a NxBJ panel of C. */
void dgepb(int N, int bj, int bk, 
	   int Astride, double A[N][Astride],  
	   int Bstride, double B[N][Bstride], 
	   int Cstride, double C[N][Cstride])
{
	int i;

	for (i = 0; i < N; i += BI) {
		int bi = min(BI, N - i);
		dgebb(bi, bj, bk, 
		      Astride, (double (*)[])&A[i][0],
		      Bstride, B,
		      Cstride, (double (*)[])&C[i][0]);
	}
}

/* Breaks dgemm down into panels of A and C and blocks of B */
void square_dgemm(int N, double A[N][N], double B[N][N], double C[N][N])
{
	int j, k;
	for (k = 0; k < N; k += BK) {
		int bk = min(BK, N - k);
		for (j = 0; j < N; j += BJ)  {
			int bj = min(BJ, N - j);
			dgepb(N, bj, bk,
			      N, (double (*)[N])&A[0][k], 
			      N, (double (*)[N])&B[k][j],
			      N, (double (*)[N])&C[0][j]);
		}
	}
}
