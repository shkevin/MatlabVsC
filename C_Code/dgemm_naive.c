#include "mex.h"
#include <math.h>
#include "matrix.h"
#include "dgemm.h"
#include <stdio.h>

/* Input Arguments */
#define A_IN    prhs[0]
#define B_IN    prhs[1]

/* Output Arguments */
#define C_OUT  plhs[0]

const char * dgemm_desc = "DGEMM with unoptimized loop index order";

void square_dgemm(int N, double A[N][N], double B[N][N], double C[N][N])
{

	int i, j, k;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			for (k = 0; k < N; k++) {
				(*C)[i + j*N] += (*A)[i+k*N] * (*B)[k+j*N];
			}
		}
	}
}


void mexFunction(int nlhs, mxArray *plhs[], 
          int nrhs, const mxArray *prhs[])
     
{ 
    double *A, *B;
    double *C;
    size_t m,n;

    // if (nrhs != 2) { 
    //     mexErrMsgIdAndTxt( "MATLAB:yprime:invalidNumInputs",
    //             "Two input arguments required."); 
    // } else if (nlhs > 1) {
    //     mexErrMsgIdAndTxt( "MATLAB:dgemm_naive:maxlhs",
    //             "Too many output arguments."); 
    // }

    // /* check to make sure the first input argument is a real matrix */
    // if( !mxIsDouble(A_IN) || mxIsComplex(A_IN)) {
    //   mexErrMsgIdAndTxt( "MATLAB:dgemm_naive:invalidA",
    //           "First input argument must be a real matrix.");
    // }

    // /* check to make sure the second input argument is a real matrix */
    // if( !mxIsDouble(B_IN) || mxIsComplex(B_IN)) {
    //   mexErrMsgIdAndTxt( "MATLAB:dgemm_naive:invalidB",
    //           "Second input argument must be a real matrix.");
    // }

    m = mxGetM(A_IN);
    n = mxGetN(A_IN);

    A = (double *)mxGetData(A_IN);
    B = (double *)mxGetData(B_IN);

    /* Create a matrix for the return argument */ 
    C_OUT = mxCreateDoubleMatrix((mwSize)n, (mwSize)n, mxREAL);

    /* Assign pointers to the various parameters */ 
    C = mxGetData(C_OUT);

    /* Do the actual computations in a subroutine */
    square_dgemm((mwSize)n,(double (*)[])A,(double (*)[])B,(double (*)[])C);
    return;
}
