//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef LAPACK_WRAPPER_HPP__
#define LAPACK_WRAPPER_HPP__


/*		Since the library is still being developed, only used functions are wrapped.
 */
namespace COPT
{

int block_size( char* name, char* opts,int n1,int n2, int n3,int n4)
{
	int ispec = 1;
	return ilaenv_(&ispec,name,opts,&n1,&n2,&n3,&n4);
}

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

template<class index,class scalar>
int copt_lapack_gelss(index m, index n, index nrhs,
	scalar* a, index lda, scalar* b,
	index ldb, scalar* s, scalar rcond,
	index* rank,
	index* info)
{
	throw COException("Unknown type for lapack wrapper!");
}
template<>
int copt_lapack_gelss(int m, int n, int nrhs,
	float *a, int lda, float *b,
	int ldb, float *s, float rcond,
	int* rank,
	int* info)
{
	int lwork = 3*std::min(m,n)+std::max(2*std::min(m,n),std::max(std::max(m,n),nrhs));
	float* work = new float[lwork];
	int status = sgelss_(&m,&n,&nrhs,a,&lda,b,&ldb,s,&rcond,rank,work,&lwork,info);
	delete[] work;
	return status;
}

template<>
int copt_lapack_gelss(int m,int n,int nrhs,
	double *a, int lda, double *b,
	int ldb, double* s,double rcond,
	int* rank,
	int* info)
{
	int lwork = 3*std::min(m,n)+std::max(2*std::min(m,n),std::max(std::max(m,n),nrhs));
	double* work = new double[lwork];
	int status = dgelss_(&m,&n,&nrhs,a,&lda,b,&ldb,s,&rcond,rank,work,&lwork,info);
	delete[] work;
	return status;
}

template<class index,class scalar>
int copt_lapack_gels(char trans, index m, index n,
	index nrhs, scalar *a, index lda,
	scalar* b, index ldb,
	index *info)
{
	throw COException("Unknow type for lapack wrapper!");
}

template<>
int copt_lapack_gels(char trans,int m, int n,
	int nrhs, float *a, int lda,
	float *b, int ldb,
	int *info)
{
	int nb = block_size("sgels","U",m,n,nrhs,-1);
	int lwork = std::min(m,n)+std::max(std::max(m,n),nrhs)*nb;
	float* work = new float[lwork];
	int status = sgels_(&trans,&m,&n,&nrhs,a,&lda,b,&ldb,work,&lwork,info);
	delete []work;
	return status;
}

template<>
int copt_lapack_gels(char trans, int m, int n,
	int nrhs, double *a, int lda,
	double *b, int ldb,
	int *info)
{
	int nb = block_size("dgels","U",m,n,nrhs,-1);
	int lwork = std::min(m,n)+std::max(std::max(m,n),nrhs)*nb;
	double* work = new double[lwork];
	int status = dgels_(&trans,&m,&n,&nrhs,a,&lda,b,&ldb,work,&lwork,info);
	delete []work;
	return status;
}


}

#endif