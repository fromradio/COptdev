//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef LAPACK_WRAPPER_HPP__
#define LAPACK_WRAPPER_HPP__


/*		Since the library is still being developed, only used functions are wrapped.
 */
namespace COPT
{

int block_size( const char* name,const char* opts,int n1,int n2, int n3,int n4)
{
	int ispec = 1;
	return ilaenv_(&ispec,const_cast<char*>(name),const_cast<char*>(opts),&n1,&n2,&n3,&n4);
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

/** symmetric matrix to tridiagonal matrix */
template<class index ,class scalar>
int copt_lapack_sytrd(char uplo , index n , scalar* a , index lda , scalar *d , scalar *e , scalar *tau ,  index *info )
{
	throw COException("Unknown type for lapack wrapper!");
}

template<class index, class scalar>
int copt_lapack_sytrd(char uplo, index n, std::complex<scalar>* a, int lda,
					scalar *d, scalar *e, std::complex<scalar>* tau, 
					index *info )
{
	throw COException("Unknown type for lapack wrapper!");
}

template<>
int copt_lapack_sytrd(char uplo , int n , float* a , int lda , float *d ,float *e , float *tau ,  int *info )
{
	/** compute the lwork first */
	float *work = new float[1];
	int lwork = -1;
	ssytrd_(&uplo,&n,a,&lda,d,e,tau,work,&lwork,info);
	lwork = work[0];
	delete[]work;
	work = new float[lwork];
	int status = ssytrd_(&uplo,&n,a,&lda,d,e,tau,work,&lwork,info);
	delete[] work;
	return status;
}

template<>
int copt_lapack_sytrd(char uplo , int n , double *a , int lda , double *d , double *e , double *tau , int *info )
{
	/** compute the lwork first */
	double *work = new double[1];
	int lwork = -1;
	dsytrd_(&uplo,&n,a,&lda,d,e,tau,work,&lwork,info);
	lwork = work[0];
	delete[]work;
	work = new double[lwork];
	int status = dsytrd_(&uplo,&n,a,&lda,d,e,tau,work,&lwork,info);
	delete[]work;
	return status;
}

template<class index,class scalar>
int copt_lapack_hetrd(char uplo, int n, scalar *a, int lda,
				scalar *d, scalar *e, scalar *tau,
				index *info )
{
	throw COException("Unknown type for lapack wrapper!");
}

template<class index,class scalar>
int copt_lapack_hetrd(char uplo, int n, std::complex<scalar>* a, int lda, 
				scalar *d, scalar *e, std::complex<scalar>* tau,
				index *info )
{
	throw COException("Unknown type for lapack wrapper!");
}

template<>
int copt_lapack_hetrd(char uplo, int n, std::complex<float> *a, int lda,
				float *d, float *e, std::complex<float> *tau,
				int *info )
{
	std::complex<float> *work = new std::complex<float>[1];
	int lwork = -1;
	chetrd_(&uplo,&n,a,&lda,d,e,tau,work,&lwork,info);
	lwork = static_cast<int>(work[0].real());
	delete[]work;
	work = new std::complex<float>[lwork];
	int status = chetrd_(&uplo,&n,a,&lda,d,e,tau,work,&lwork,info);
	delete[]work;
	return status;
}

template<>
int copt_lapack_hetrd(char uplo, int n, std::complex<double> *a, int lda,
				double *d, double *e, std::complex<double> *tau,
				int *info )
{
	std::complex<double> *work = new std::complex<double>[1];
	int lwork = -1;
	zhetrd_(&uplo,&n,a,&lda,d,e,tau,work,&lwork,info);
	lwork = static_cast<int>(work[0].real());
	delete[]work;
	work = new std::complex<double>[lwork];
	int status = zhetrd_(&uplo,&n,a,&lda,d,e,tau,work,&lwork,info);
	delete[] work;
	return status;
}

/** stebz */
template<class index,class scalar>
int copt_lapack_stebz(char range,char order,index n,scalar vl,scalar vu,index il,index iu,scalar abstol, scalar* d, scalar *e , index *m , index *nsplit , scalar *w , index *iblock , index *isplit , scalar *work , index *iwork , index *info )
{
	throw COException("Unknown type for lapack wrapper!");
}

template<>
int copt_lapack_stebz(char range,char order,int n , float vl , float vu , int il , int iu , float abstol , float *d , float * e , int *m , int* nsplit , float *w , int* iblock , int* isplit , float* work, int* iwork , int* info)
{
	return sstebz_(&range,&order,&n,&vl,&vu,&il,&iu,&abstol,d,e,m,nsplit,w,iblock,isplit,work,iwork,info);
}

template<>
int copt_lapack_stebz(char range , char order , int n , double vl , double vu , int il , int iu , double abstol , double *d , double *e , int *m , int *nsplit , double *w , int *iblock , int *isplit , double *work , int *iwork , int *info )
{
	return dstebz_(&range,&order,&n,&vl,&vu,&il,&iu,&abstol,d,e,m,nsplit,w,iblock,isplit,work,iwork,info);
}

/** LU factorization */
template<class index,class scalar>
int copt_lapack_getrf( index m , index n , scalar *a,
				index lda , index *ipiv,
				index *info )
{
	throw COException("Unkonwn type for lapack wrapper!");
}

template<>
int copt_lapack_getrf( int m , int n , float *a,
				int lda , int *ipiv ,
				int *info )
{
	return sgetrf_(&m,&n,a,&lda,ipiv,info);
}

template<>
int copt_lapack_getrf( int m , int n , double *a,
				int lda , int *ipiv ,
				int *info )
{
	return dgetrf_(&m,&n,a,&lda,ipiv,info);
}

template<>
int copt_lapack_getrf( int m, int n, std::complex<float> *a,
				int lda, int *ipiv,
				int *info )
{
	return cgetrf_(&m,&n,a,&lda,ipiv,info);
}

template<>
int copt_lapack_getrf( int m, int n, std::complex<double> *a,
				int lda, int *ipiv,
				int *info)
{
	return zgetrf_(&m,&n,a,&lda,ipiv,info);
}

/** LU solving */
template<class index,class scalar>
int copt_lapack_getrs( char trans , index n , index nrhs,
				scalar *a, index lda , index * ipiv,
				scalar *b, index ldb,
				index *info)
{
	throw COException("Unkonw type for lapack wrapper!");
}

template<>
int copt_lapack_getrs( char trans , int n , int nrhs,
				float *a , int lda , int *ipiv,
				float *b , int ldb ,
				int *info)
{
	return sgetrs_(&trans,&n,&nrhs,a,&lda,ipiv,b,&ldb,info);
}

template<>
int copt_lapack_getrs( char trans , int n , int nrhs,
				double *a , int lda , int *ipiv,
				double *b , int ldb ,
				int *info)
{
	return dgetrs_(&trans,&n,&nrhs,a,&lda,ipiv,b,&ldb,info);
}

template<>
int copt_lapack_getrs( char trans, int n, int nrhs,
				std::complex<float> *a, int lda, int *ipiv,
				std::complex<float> *b, int ldb,
				int *info )
{
	return cgetrs_(&trans,&n,&nrhs,a,&lda,ipiv,b,&ldb,info);
}

template<>
int copt_lapack_getrs( char trans, int n, int nrhs,
				std::complex<double> *a, int lda, int *ipiv,
				std::complex<double> *b, int ldb,
				int *info )
{
	return zgetrs_(&trans,&n,&nrhs,a,&lda,ipiv,b,&ldb,info);
}

/** LU inverse */
template<class index,class scalar>
int copt_lapack_getri( index n , scalar *a , index lda ,
				index *ipiv , index *info )
{
	throw COException("Unkown type for lapack wrapper!");
}

template<>
int copt_lapack_getri( int n , float *a , int lda ,
				int *ipiv , int *info )
{
	float *work = new float[1];
	int lwork = -1;
	sgetri_(&n,a,&lda,ipiv,work,&lwork,info);
	lwork = work[0];
	delete[]work;
	work = new float[lwork];
	int status = sgetri_(&n,a,&lda,ipiv,work,&lwork,info);
	delete[]work;
	return status;
}

template<>
int copt_lapack_getri( int n , double *a , int lda , 
				int *ipiv , int *info )
{
	double *work = new double[1];
	int lwork = -1;
	dgetri_(&n,a,&lda,ipiv,work,&lwork,info);
	lwork = work[0];
	delete[] work;
	work = new double[lwork];
	int status = dgetri_(&n,a,&lda,ipiv,work,&lwork,info);
	delete[]work;
	return status;
}

template<>
int copt_lapack_getri( int n, std::complex<float> *a, int lda,
				int *ipiv, int *info )
{
	std::complex<float> *work = new std::complex<float>[1];
	int lwork = -1;
	cgetri_(&n,a,&lda,ipiv,work,&lwork,info);
	lwork = static_cast<int>(work[0].real());
	delete[] work;
	work = new std::complex<float>[lwork];
	int status = cgetri_(&n,a,&lda,ipiv,work,&lwork,info);
	delete[]work;
	return status;
}

template<>
int copt_lapack_getri( int n, std::complex<double> *a ,int lda,
				int *ipiv, int *info )
{
	std::complex<double> *work = new std::complex<double>[1];
	int lwork = -1;
	zgetri_(&n,a,&lda,ipiv,work,&lwork,info);
	lwork = static_cast<int>(work[0].real());
	delete[] work;
	work = new std::complex<double>[lwork];
	int status = zgetri_(&n,a,&lda,ipiv,work,&lwork,info);
	delete[]work;
	return status;
}

/** Cholesky factorization */
template<class index,class scalar>
int copt_lapack_potrf( char uplo , index n , scalar * a,
				index lda ,
				index *info )
{
	throw COException("Unknown type for lapack wrapper!");
}

template<>
int copt_lapack_potrf( char uplo , int n , float *a ,
				int lda,
				int *info )
{
	return spotrf_(&uplo,&n,a,&lda,info);
}

template<>
int copt_lapack_potrf( char uplo , int n , double *a ,
				int lda,
				int *info )
{
	return dpotrf_(&uplo,&n,a,&lda,info);
}

template<>
int copt_lapack_potrf( char uplo, int n, std::complex<float> *a,
				int lda,
				int *info)
{
	return cpotrf_(&uplo,&n,a,&lda,info);
}

template<>
int copt_lapack_potrf( char uplo, int n, std::complex<double> *a,
				int lda,
				int *info)
{
	return zpotrf_(&uplo,&n,a,&lda,info);
}

/** Cholesky solving */
template<class index,class scalar>
int copt_lapack_potrs(char uplo , index n , index nrhs,
				scalar *a , index lda , scalar *b,
				index ldb,
				index *info )
{
	throw COException("Unknown type for lapack wrapper!");
}

template<>
int copt_lapack_potrs(char uplo , int n , int nrhs,
				float *a , int lda , float *b,
				int ldb,
				int *info )
{
	return spotrs_(&uplo,&n,&nrhs,a,&lda,b,&ldb,info);
}

template<>
int copt_lapack_potrs(char uplo , int n , int nrhs,
				double *a , int lda , double *b,
				int ldb,
				int *info)
{
	return dpotrs_(&uplo,&n,&nrhs,a,&lda,b,&ldb,info);
}

template<>
int copt_lapack_potrs( char uplo, int n, int nrhs,
				std::complex<float> *a, int lda, std::complex<float> *b,
				int ldb,
				int *info)
{
	return cpotrs_(&uplo,&n,&nrhs,a,&lda,b,&ldb,info);
}

template<>
int copt_lapack_potrs( char uplo, int n, int nrhs,
				std::complex<double> *a, int lda, std::complex<double> *b,
				int ldb,
				int *info )
{
	return zpotrs_(&uplo,&n,&nrhs,a,&lda,b,&ldb,info);
}

/** Cholesky inverse */
template<class index,class scalar>
int copt_lapack_potri( char uplo , index n , scalar *a,
				int lda,
				int *info)
{
	throw COException("Unknown type for lapack wrapper!");
}

template<>
int copt_lapack_potri( char uplo , int n , float *a,
				int lda,
				int *info)
{
	return spotri_(&uplo,&n,a,&lda,info);
}

template<>
int copt_lapack_potri( char uplo , int n , double *a,
				int lda,
				int *info)
{
	return dpotri_(&uplo,&n,a,&lda,info);
}

template<>
int copt_lapack_potri( char uplo, int n, std::complex<float> *a,
				int lda,
				int *info )
{
	return cpotri_(&uplo,&n,a,&lda,info);
}

template<>
int copt_lapack_potri( char uplo, int n, std::complex<double> *a,
				int lda,
				int *info )
{
	return zpotri_(&uplo,&n,a,&lda,info);
}

/** qr factorization */
template<class index,class scalar>
int copt_lapack_geqrf( index m , index n , scalar *a,
				index lda , scalar *tau,
				index *info)
{
	throw COException("Unknown type for lapack wrapper!");
}

template<>
int copt_lapack_geqrf( int m , int n , float *a,
				int lda , float *tau,
				int *info)
{
	float *work = new float[1];
	int lwork = -1;
	sgeqrf_(&m,&n,a,&lda,tau,work,&lwork,info);
	lwork = work[0];
	delete[] work;
	work = new float[lwork];
	int status = sgeqrf_(&m,&n,a,&lda,tau,work,&lwork,info);
	delete[]work;
	return status;
}

template<>
int copt_lapack_geqrf( int m , int n , double *a,
				int lda , double *tau ,
				int *info )
{
	double *work = new double[1];
	int lwork = -1;
	dgeqrf_(&m,&n,a,&lda,tau,work,&lwork,info);
	lwork = work[0];
	delete[] work;
	work = new double[lwork];
	int status = dgeqrf_(&m,&n,a,&lda,tau,work,&lwork,info);
	delete[]work;
	return status;
}

template<>
int copt_lapack_geqrf( int m , int n , std::complex<float> *a,
				int lda , std::complex<float> *tau,
				int *info )
{
	std::complex<float> *work = new std::complex<float>[1];
	int lwork = -1;
	cgeqrf_(&m,&n,a,&lda,tau,work,&lwork,info);
	lwork = static_cast<int>(work[0].real());
	delete[] work;
	work = new std::complex<float>[lwork];
	int status = cgeqrf_(&m,&n,a,&lda,tau,work,&lwork,info);
	delete[] work;
	return status;
}

template<>
int copt_lapack_geqrf( int m , int n , std::complex<double> *a,
				int lda , std::complex<double>* tau,
				int *info)
{
	std::complex<double> *work = new std::complex<double>[1];
	int lwork = -1;
	zgeqrf_(&m,&n,a,&lda,tau,work,&lwork,info);
	lwork = static_cast<int>(work[0].real());
	delete[] work;
	work = new std::complex<double>[lwork];
	int status = zgeqrf_(&m,&n,a,&lda,tau,work,&lwork,info);
	delete[]work;
	return status;
}

/** compute Q'*B */
template<class index , class scalar>
int copt_lapack_ormqr(char side,char trans,index m , index n,
				index k,scalar *a , index lda , scalar *tau , scalar *c , index ldc,
				index *info)
{
	throw COException("Unknown type for lapack wrapper!");
}

template<>
int copt_lapack_ormqr(char side,char trans,int m,int n,
				int k,float *a,int lda,float *tau,float *c,int ldc,
				int *info)
{
	float *work = new float[1];
	int lwork = -1;
	sormqr_(&side,&trans,&m,&n,&k,a,&lda,tau,c,&ldc,work,&lwork,info);
	lwork = work[0];
	delete[] work;
	work = new float[lwork];
	int status = sormqr_(&side,&trans,&m,&n,&k,a,&lda,tau,c,&ldc,work,&lwork,info);
	delete[] work;
	return status;
}

template<>
int copt_lapack_ormqr(char side,char trans,int m,int n,
				int k, double *a, int lda , double *tau, double *c, int ldc,
				int *info)
{
	double *work = new double[1];
	int lwork = -1;
	dormqr_(&side,&trans,&m,&n,&k,a,&lda,tau,c,&ldc,work,&lwork,info);
	lwork = work[0];
	delete[] work;
	work = new double[lwork];
	int status = dormqr_(&side,&trans,&m,&n,&k,a,&lda,tau,c,&ldc,work,&lwork,info);
	delete[] work;
	return status;
}

template<>
int copt_lapack_ormqr( char side,char trans, int m, int n,
				int k, std::complex<float> *a, int lda, std::complex<float> *tau,
				std::complex<float> *c, int ldc,
				int *info)
{
	std::complex<float> *work = new std::complex<float>[1];
	int lwork = -1;
	cunmqr_(&side,&trans,&m,&n,&k,a,&lda,tau,c,&ldc,work,&lwork,info);
	lwork = static_cast<int>(work[0].real());
	delete[] work;
	work = new std::complex<float>[lwork];
	int status = cunmqr_(&side,&trans,&m,&n,&k,a,&lda,tau,c,&ldc,work,&lwork,info);
	return status;
}

template<>
int copt_lapack_ormqr( char side, char trans, int m, int n,
				int k, std::complex<double> *a, int lda, std::complex<double> *tau,
				std::complex<double> *c, int ldc,
				int *info )
{
	std::complex<double> *work = new std::complex<double>[1];
	int lwork = -1;
	zunmqr_(&side,&trans,&m,&n,&k,a,&lda,tau,c,&ldc,work,&lwork,info);
	lwork = static_cast<int>(work[0].real());
	delete[] work;
	work = new std::complex<double>[lwork];
	int status = zunmqr_(&side,&trans,&m,&n,&k,a,&lda,tau,c,&ldc,work,&lwork,info);
	return status;
}

/** strsm */

// template<class index,class scalar>
// int copt_lapack_trsm(char side,char uplo,char transa,char diag,
// 				index m , index n , scalar alpha, scalar *a, index lda , scalar *b,
// 				index ldb )
// {
// 	throw COException("Unknown type for lapack wrapper!");
// }

// template<>
int copt_lapack_trsm(char side,char uplo, char transa, char diag,
				int m, int n, float alpha, float *a, int lda, float *b,
				int ldb )
{
	return strsm_(&side,&uplo,&transa,&diag,&m,&n,&alpha,a,&lda,b,&ldb);
}

// template<>
int copt_lapack_trsm(char side,char uplo,char transa,char diag,
				int m, int n, double alpha, double *a, int lda, double *b,
				int ldb )
{
	return dtrsm_(&side,&uplo,&transa,&diag,&m,&n,&alpha,a,&lda,b,&ldb);
}

int copt_lapack_trsm(char side,char uplo, char transa, char diag,
				int m, int n, std::complex<float> alpha, std::complex<float> *a, int lda,
				std::complex<float> *b, int ldb )
{
	return ctrsm_(&side,&uplo,&transa,&diag,&m,&n,&alpha,a,&lda,b,&ldb);
}

int copt_lapack_trsm(char side, char uplo, char transa, char diag,
				int m, int n, std::complex<double> alpha, std::complex<double> *a, int lda,
				std::complex<double> *b, int ldb )
{
	return ztrsm_(&side,&uplo,&transa,&diag,&m,&n,&alpha,a,&lda,b,&ldb);
}

/** solve QR problem */
template<class index,class scalar>
int copt_lapack_geqrs( index m , index n , index nrhs , 
				scalar *a , index lda , scalar *tau, scalar *b, index ldb,
				index *info)
{
	throw COException("Unknown type for lapack wrapper!");
}

template<>
int copt_lapack_geqrs(int m , int n , int nrhs ,
				float *a , int lda , float *tau , float *b , int ldb,
				int *info )
{
	// compute B: = Q'*B
	copt_lapack_ormqr('L','T',m,nrhs,n,a,lda,tau,b,ldb,info);
	// solve R*x = B
	return copt_lapack_trsm('L','U','N','N',n,nrhs,1.0,a,lda,b,ldb);
}

template<>
int copt_lapack_geqrs(int m, int n, int nrhs,
				double *a, int lda , double *tau, double *b , int ldb,
				int *info)
{
	// compute B:= Q'*B
	copt_lapack_ormqr('L','T',m,nrhs,n,a,lda,tau,b,ldb,info);
	// solve R*x = B
	return copt_lapack_trsm('L','U','N','N',n,nrhs,1.0,a,lda,b,ldb);
}

template<>
int copt_lapack_geqrs(int m, int n, int nrhs,
				std::complex<float> *a, int lda, std::complex<float> *tau,
				std::complex<float> *b, int ldb,
				int *info )
{
	// compute B:= Q'*B
	copt_lapack_ormqr('L','C',m,nrhs,n,a,lda,tau,b,ldb,info);
	// solve R*x = B
	return copt_lapack_trsm('L','U','N','N',n,nrhs,std::complex<float>(1.0,0),a,lda,b,ldb);
}

template<>
int copt_lapack_geqrs(int m, int n,int nrhs,
				std::complex<double> *a, int lda, std::complex<double> *tau,
				std::complex<double> *b, int ldb,
				int *info )
{
	// compute B:= Q'*B
	copt_lapack_ormqr('L','C',m,nrhs,n,a,lda,tau,b,ldb,info);
	// solve R*x = B
	return copt_lapack_trsm('L','U','N','N',n,nrhs,std::complex<double>(1.0,0),a,lda,b,ldb);
}

}

#endif