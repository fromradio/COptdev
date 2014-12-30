#include "../include/Header"
#include "mex.h"

#define A prhs[0]
#define b prhs[1]

typedef COPT::KernelTrait<double,int> 		kernel;
typedef kernel::Matrix						Matrix;
typedef kernel::Vector 						Vector;
void mexFunction(int nlhs, mxArray *plhs[],
				int nrhs, const mxArray *prhs[])
{
	int M,N,NN,K;
	double *a_arr,*b_arr;
	if ( nrhs != 2 )
		mexErrMsgTxt("Wrong number of input arguments!");
	else if ( nlhs > 1 )
		mexErrMsgTxt("Too many output arguments!");

	M = mxGetM(A);
	N = mxGetN(A);
	NN = mxGetM(b);
	K = mxGetN(b);
	a_arr = mxGetPr(A);
	b_arr = mxGetPr(b);

	if (M!=NN)
		mexErrMsgTxt("Wrong input!");
	Matrix AL(M,N,a_arr);
	Matrix br(NN,K,b_arr);

	COPT::LU<Matrix> lu(AL);
	Matrix x=lu.solve(br);

	plhs[0] = mxCreateDoubleMatrix(N,K,mxREAL);
	double *l = mxGetPr(plhs[0]);
	for ( int i = 0 ; i < x.size() ; ++ i )
	{
		l[i] = x.dataPtr()[i];
	}

	return;
}