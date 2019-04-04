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

#if (VECTLEN == 2)
#define VPREFIX _mm128
#define VTYPE __m128d
#define VOP(op) _mm_ ## op
#define VBROADCAST _mm_loaddup_pd
#define VSHIFT 1
#elif (VECTLEN == 4)
#define VPREFIX _mm256
#define VTYPE __m256d
#define VOP(op) _mm256_
#define VBROADCAST __mm256_broadcast_sd
#define VSHIFT 2
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
	VTYPE a, blocal[RJ>>VSHIFT], clocal[RI][RJ>>VSHIFT];
	double (*Bin)[Bstride] __attribute__ ((aligned (VECTLEN<<3)));
	int i, j, k;

	Bin = __builtin_assume_aligned(B, VECTLEN<<3);
        for (i = 0; i < RI; i++)
                for (j = 0; j < (RJ>>VSHIFT); j++)
                        clocal[i][j] = VOP(loadu_pd)(&C[i][j<<VSHIFT]);

        for (k = 0; k < bk; k++) {
                for (j = 0; j < (RJ>>VSHIFT); j++) {
                        blocal[j] = VOP(load_pd)(&Bin[k][j<<VSHIFT]);
                }
                for (i = 0; i < RI; i++) {
                        a = VBROADCAST(&A[i][k]);
                        for (j = 0; j < (RJ>>VSHIFT); j++) {
                                clocal[i][j] = VOP(add_pd)(clocal[i][j],
                                                             VOP(mul_pd)(a, blocal[j]));
                        }
                }
        }
        for (i = 0; i < RI; i++) {
                for (j = 0; j < (RJ>>VSHIFT); j++) {
                        VOP(storeu_pd)(&C[i][j<<VSHIFT], clocal[i][j]);
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
