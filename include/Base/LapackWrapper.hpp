//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef LAPACK_WRAPPER_HPP__
#define LAPACK_WRAPPER_HPP__


/*		Since the library is still being developed, only used functions are wrapped.
 */
namespace COPT
{

template<class index,class real>
int copt_lapack_gesv( index *n, index *nrhs, real *a,
	index *lda, index *ipiv, real *b,
	index *ldb,
	index *info)
{
	throw COException("unknown type for lapack wrapper!");
	return -1;
}
template<>
int copt_lapack_gesv(int *n, int *nrhs, float *a,
	int *lda, int *ipiv, float *b,
	int *ldb,
	int* info)
{
	return sgesv_(n,nrhs,a,lda,ipiv,b,ldb,info);
}

template<>
int copt_lapack_gesv(int *n, int *nrhs,	 double *a,
	int *lda, int *ipiv, double *b,
	int *ldb,
	int *info)
{
	return dgesv_(n,nrhs,a,lda,ipiv,b,ldb,info);
}

template<>
int copt_lapack_gelss(int m, int n, int nrhs,
	float *a, int lda, float *b,
	int ldb, float *s, float *rcond,
	int* rank, float *work, int *lwork,
	int* info)
{
	return sgelss_(&m,&n,&nrhs,a,&lda,b,&ldb,s,rcond,rank,work,lwork,info);
}

template<>
int copt_lapack_gelss(int m, int n, int nrhs,
	double *a, int lda, double *b,
	int ldb, double *s, double *rcond,
	int *rank, double *work, int *lwork,
	int *info)
{
	return dgelss_(&m,&n,&nrhs,a,&lda,b,&ldb,s,rcond,rank,work,lwork,info);
}


}

#endif