#include "dgemm.h"
#include <stdio.h>

#include <immintrin.h>
#include <emmintrin.h>

#define xstr(s) str(s)
#define str(s) #s

static inline int min(int i, int j)
{
	return i < j ? i : j;
}


/* Register blocking parameters */
#ifndef RI
#define RI 2
#endif

#ifndef RJ
#define RJ 16
#endif

void dgebb_subblock_opt(int bk,
			int Astride, double A[][Astride],  
	  	        int Bstride, double B[][Bstride], 
	 	        int Cstride, double C[][Cstride])
{
	__m256d alocal[RI], blocal[RJ>>2], clocal[RI][RJ>>2];
	int i, j, k;
	for (i = 0; i < RI; i++) {
		for (j = 0; j < (RJ>>2); j++) {
			clocal[i][j] = _mm256_loadu_pd(&C[i][j<<2]);
		}
	}

	for (k = 0; k < bk; k++) {
		for (i = 0; i < RI; i++) {
			alocal[i] = _mm256_broadcast_sd(&A[i][k]);
		}
		for (j = 0; j < (RJ>>2); j++) {
 			blocal[j] = _mm256_load_pd(&B[k][j<<2]);
		}
		for (i = 0; i < RI; i++) {
			for (j = 0; j < (RJ>>2); j++) {
				__m256d tmp = _mm256_mul_pd(alocal[i], 
							    blocal[j]);
				clocal[i][j] = _mm256_add_pd(clocal[i][j], tmp);
							     
			}
		}
	}
	for (i = 0; i < RI; i++) {
		for (j = 0; j < (RJ>>2); j++) {
			_mm256_storeu_pd(&C[i][j<<2], clocal[i][j]);
		}
	}
}

void dgebb_subblock_gen(int bi, int bj, int bk,
			int Astride, double A[][Astride],  
			int Bstride, double B[][Bstride], 
			int Cstride, double C[][Cstride])
{
	int i, k, j;
	for (i = 0; i < bi; i++) {
		for (k = 0; k < bk; k++) {
			double tmp = A[i][k];
			for (j = 0; j < bj; j++) {
				C[i][j] += tmp * B[k][j];
			}
		}
	}
}


/* Multiply a bi*bk block of A times a bk*bj block of B 
 * into an bi * bj block of C. Subblock for registers when possible */
void dgebb(int bi, int bj, int bk,
	   int Astride, double A[][Astride],  
	   int Bstride, double B[][Bstride], 
	   int Cstride, double C[][Cstride])
{
	/*
	 * This code goes through all blocks that can be optimized first, 
	 * then cleans up the edges. Removes some conditional overheads
	 * at the outer level for modest performance gains */
	int i, j;
	int iend = 0, jend = 0;

	/* Get the main body of the block */
	for (i = 0; i < bi + 1 - RI; i+=RI) {
		for (j = 0; j < bj + 1 - RJ; j+= RJ) {
			dgebb_subblock_opt(bk, 
					   Astride, (double (*)[])&A[i][0],
					   Bstride, (double (*)[])&B[0][j],
		       		   	   Cstride, (double (*)[])&C[i][j]);
		}
	}
	iend = bi - bi % RI;
	jend = bj - bj % RJ;

	/* There are some rows we didn't cover */
	if (iend < bi) {
		dgebb_subblock_gen(bi - iend, jend, bk,
			           Astride, (double (*)[])&A[iend][0],
				   Bstride, (double (*)[])&B[0][0],
			           Cstride, (double (*)[])&C[iend][0]);
	} 
	/* There are some columns we didn't cover */
	if (jend < bj) {
		dgebb_subblock_gen(iend, bj - jend, bk,
			           Astride, (double (*)[])&A[0][0],
				   Bstride, (double (*)[])&B[0][jend],
			           Cstride, (double (*)[])&C[0][jend]);
	} 
	/* There are some row/column pairs we didn't cover */
	if (iend < bi && jend < bj) {
		dgebb_subblock_gen(bi - iend, bj - jend, bk,
			           Astride, (double (*)[])&A[iend][0],
				   Bstride, (double (*)[])&B[0][jend],
			           Cstride, (double (*)[])&C[iend][jend]);
	}
}


/* Cache blocking parameters */
#ifndef BI
#define BI 128
#endif

#ifndef BJ
#define BJ 32
#endif

#ifndef BK
#define BK 64
#endif

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

#define CBLOCK_DESCR BI * BJ * BK
#define RBLOCK_DESCR RI * RJ

const char * dgemm_desc = "Vectorized DGEMM packed/cache blocked " xstr(CBLOCK_DESCR) 
			   "/register blocked " xstr(RBLOCK_DESCR) " with optimized DGEBB";

/* Breaks dgemm down into panels of A and C and blocks of B */
void square_dgemm(int N, double A[N][N], double B[N][N], double C[N][N])
{
	int j, k;
	for (k = 0; k < N; k += BK) {
		int bk = min(BK, N - k);
		for (j = 0; j < N; j += BJ)  {
			int bj = min(BJ, N - j);
			int jj, kk;
			/* Copy the key block into contiguous local memory */
			static double blocal[BK][BJ]  
			       __attribute__ ((aligned (32)));
			for (kk = 0; kk < bk; kk++) {
				for (jj = 0; jj < bj; jj++) {
					blocal[kk][jj] = B[k+kk][j+jj];
				}
			}
			dgepb(N, bj, bk,
			      N, (double (*)[N])&A[0][k], 
			      BJ, (double (*)[N])&blocal[0][0],
			      N, (double (*)[N])&C[0][j]);
		}
	}
}
