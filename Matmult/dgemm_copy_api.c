#include "dgemm.h"
#include "membuffer.h"

#include <stdio.h>

#include <immintrin.h>
#include <emmintrin.h>

#define PACKB

/* Cache blocking parameters */
#ifndef BI
#define BI 32
#endif

#ifndef BJ
#define BJ 32
#endif

#ifndef BK
#define BK 32
#endif

/* Register blocking parameters */
#ifndef RI
#define RI 2
#endif

#ifndef RJ
#define RJ 2
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
	double a, blocal[RJ], clocal[RI][RJ];
	int i, j, k;
	for (i = 0; i < RI; i++)
		for (j = 0; j < RJ; j++)
			clocal[i][j] = C[i][j];

	for (k = 0; k < bk; k++) {
		for (j = 0; j < RJ; j++) {
 			blocal[j] = B[k][j];
		}
		for (i = 0; i < RI; i++) {
			a = A[i][k];
			for (j = 0; j < RJ; j++) {
				clocal[i][j] = clocal[i][j] + a * blocal[j];
			}
		}
	}
	for (i = 0; i < RI; i++) {
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

void copy_routine(void *dest, void *src, void *arg)
{
        int jj, kk;
	int *sizes = (int *)arg;
	int N = sizes[0];
	int bk = sizes[1], bj = sizes[2];

        /* Copy the key block into contiguous local memory */
	double (*pin)[N], (*pout)[BJ];
	pin = src; pout = dest;
        for (kk = 0; kk < bk; kk++) {
                for (jj = 0; jj < bj; jj++) {
                        pout[kk][jj] = pin[kk][jj];
                }
        }
}

/* Breaks dgemm down into panels of A and C and blocks of B */
void square_dgemm(int N, double A[N][N], double B[N][N], double C[N][N])
{
	int j, k;
	struct membuffer mbhandle;

	membuffer_init(&mbhandle, 0, BK*BJ*sizeof(double), 1, copy_routine, NULL);

	for (k = 0; k < N; k += BK) {
		int bk = min(BK, N - k);
		for (j = 0; j < N; j += BJ)  {
			int bj = min(BJ, N - j);
#ifdef PACKB
			double (*bp)[BJ];
			int sizes[3];
			sizes[0] = N; sizes[1] = bk; sizes[2] = bj;
			membuffer_start_copyin(&mbhandle, 0, &B[k][j], sizes);
			bp = (double (*)[BJ])membuffer_get_block(&mbhandle, 0);
#endif

			dgepb(N, bj, bk,
			      N, (double (*)[N])&A[0][k], 
#ifdef PACKB
                              BJ, (double (*)[N])bp,
#else
                              N, (double (*)[N])&B[k][j],
#endif
			      N, (double (*)[N])&C[0][j]);
			membuffer_release_block(&mbhandle, 0);
		}
	}
	membuffer_deinit(&mbhandle);
}
