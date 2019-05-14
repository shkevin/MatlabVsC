/* THis code is strutured so that a smart compiler (read the Intel compiler)
 * will vectorize it as effetively as the hand-vectorized one in genvect. Note
 * that GCC 5.4.0 doesn't successfully vectorize this code well. I haven't 
 * checked GCC 7.X.X; wit the appropriate options it may be able to do
 * what the Intel compiler does with this code. It's not structured ot 
 * ectorize with fused multiply/add since that's not available on
 * many Intel vecotr instructions sets, increasing its portability but
 * reducing its performance on the newest hardware.
 *
 * Intel compiler options needed: -O3 -qopenmp -mkl
 */

#include "dgemm.h"

#include <stdio.h>

#include <immintrin.h>
#include <emmintrin.h>

#if defined(__AVX__)
#define VECTLEN 2
#define NREGS 32
#elif defined(__SSE2__)
/* Default to SSE vectors if nothing else */
#define VECTLEN 2
#define NREGS 16
#else
#error UNKNOWN VECTOR ARCHITECTURE
#endif

/* Register blocking parameters */
/* Our main loop is register blocked for half of the registers in clocal, 
 * a quarter in blocal, and the rest for the compiler as temps. We get that
 * by making RI 2 and RJ big enough that RI*RJ is half of the registers */
#ifndef RI
#define RI 2
#endif

#ifndef RJ
#define RJ (VECTLEN*NREGS/4)
#endif

/* Cache blocking parameters */
#ifndef BI
#define BI 96
#endif

#ifndef BJ
#define BJ 32
#endif

#ifndef BK
#define BK 32
#endif

#define xstr(s) str(s)
#define str(s) #s

#define CBLOCK_DESCR BI * BJ * BK
#define RBLOCK_DESCR RI * RJ


static inline int min(int i, int j)
{
	return i < j ? i : j;
}

const char * dgemm_desc = "DGEMM cache blocked/packed " xstr(CBLOCK_DESCR) " and register blocked " xstr(RBLOCK_DESCR);

void dgebb_subblock_opt(int bk,
			int Astride, double A[][Astride],  
	  	        int Bstride, double B[][Bstride], 
	 	        int Cstride, double C[][Cstride])
{
	double a __attribute__ ((aligned (VECTLEN<<3))),
	       blocal[RJ]__attribute__ ((aligned (VECTLEN<<3))), 
	       clocal[RI][RJ] __attribute__ ((aligned (VECTLEN<<3)));
	double (*Bin)[Bstride] __attribute__ ((aligned (VECTLEN<<3)));
	int i, j, k;

	Bin = __builtin_assume_aligned(B, VECTLEN<<3);
	for (i = 0; i < RI; i++)
#pragma omp simd
		for (j = 0; j < RJ; j++)
			clocal[i][j] = C[i][j]; /* These can be unaligned */

	for (k = 0; k < bk; k++) {
#pragma omp simd
		for (j = 0; j < RJ; j++) {
 			blocal[j] = Bin[k][j];
		}
		for (i = 0; i < RI; i++) {
			a = A[i][k];
#pragma omp simd
			for (j = 0; j < RJ; j++) {
				clocal[i][j] = clocal[i][j] + a * blocal[j];
			}
		}
	}
	for (i = 0; i < RI; i++) {
#pragma omp simd
		for (j = 0; j < RJ; j++) {
			C[i][j] = clocal[i][j];
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
        int i, j;

        /* Get the main body of the block */
        for (i = 0; i < bi; i+=RI) {
                int bii = min(RI, bi - i);
                for (j = 0; j < bj; j+= RJ) {
                        int bjj = min(RJ, bj - j);
                        if ((bii == RI) && (bjj == RJ)) {
                                dgebb_subblock_opt(bk,
                                                Astride, (double (*)[])&A[i][0],
                                                Bstride, (double (*)[])&B[0][j],
                                                Cstride, (double (*)[])&C[i][j]);
                        } else {
                                dgebb_subblock_gen(bii, bjj, bk,
                                                Astride, (double (*)[])&A[i][0],
                                                Bstride, (double (*)[])&B[0][j],
                                                Cstride, (double (*)[])&C[i][j]);
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
