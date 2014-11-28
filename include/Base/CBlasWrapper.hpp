// Copyright (C) 2014 Ruimin Wang, ruimin.wang13@gmail.com
// Copyright (C) 2014 MathU



#ifndef CBLAS_WRAPPER_H
#define CBLAS_WRAPPER_H

// 		This file defines C++ wrapper for library cblas
#ifdef COPT_USE_BLAS
#endif

#ifdef CBLAS
namespace blas{
//		level 1
template<class eT>
void copt_blas_copy( const int N,const eT* X,const int incX, eT* Y,const int incY)
{
	for ( int i = 0 ; i < N ; ++ i ){
		Y[i*incY]=X[i*incX];
	}
}

void copt_blas_copy( const int N , const double* X, const int incX , double* Y , const int incY)
{
	cblas_dcopy(N,X,incX,Y,incY);
}

void copt_blas_copy( const int N, const float* X , const int incX , float* Y , const int incY )
{
	cblas_scopy(N,X,incX,Y,incY);
}

void copt_blas_copy(const int N, const std::complex<float>* X,const int incX, std::complex<float>* Y,const int incY)
{
	cblas_ccopy(N,X,incX,Y,incY);
}

void copt_blas_copy(const int N, const std::complex<double>* X,const int incX,std::complex<double>* Y,const int incY)
{
	cblas_zcopy(N,X,incX,Y,incY);
}

/*			overload of swap operation
 *
 *
 */
template<class eT>
void copt_blas_swap (const int N,eT* X,const int incX,eT* Y,const int incY)
{
}
void copt_blas_swap (const int N,float* X,const int incX,float* Y,const int incY)
{
	cblas_sswap(N,X,incX,Y,incY);
}
void copt_blas_swap (const int N,double* X,const int incX,double* Y,const int incY)
{
	cblas_dswap(N,X,incX,Y,incY);
}
void copt_blas_swap (const int N,std::complex<float>* X,const int incX,std::complex<float>* Y,const int incY)
{
	cblas_cswap(N,X,incX,Y,incY);
}
void copt_blas_swap (const int N,std::complex<double>* X,const int incX,std::complex<double>* Y,const int incY)
{
	cblas_zswap(N,X,incX,Y,incY);
}

/*
 *				dot operation
 */
template<class eT>
eT copt_blas_dot(const int N,const eT* X,const int incX,const eT* Y,const int incY)
{
	return static_cast<eT>(0.0);
}
double copt_blas_dot(const int N,const double* X,const int incX,const double* Y,const int incY)
{
	return cblas_ddot(N,X,incX,Y,incY);
}
float copt_blas_dot(const int N,const float* X,const int incX,const float* Y,const int incY)
{
	return cblas_sdot(N,X,incX,Y,incY);
}
std::complex<float> copt_blas_dot(const int N,const std::complex<float>* X,const int incX,const std::complex<float>* Y,int incY)
{
	std::complex<float> re;
	cblas_cdotu_sub(N,X,incX,Y,incY,&re);
	return re;
}
std::complex<double> copt_blas_dot(const int N,const std::complex<double>* X,const int incX,const std::complex<double>* Y,int incY)
{
	std::complex<double> re;
	cblas_zdotu_sub(N,X,incX,Y,incY,&re);
	return re;
}

/*				scale operation
 */
template<class eT>
void copt_blas_scal(const int N,const eT alpha,eT* X,const int incX)
{
}
void copt_blas_scal(const int N,const float alpha,float* X,const int incX)
{
	cblas_sscal(N,alpha,X,incX);
}
void copt_blas_scal(const int N,const double alpha,double* X,const int incX)
{
	cblas_dscal(N,alpha,X,incX);
}
void copt_blas_scal(const int N,const std::complex<float>& alpha,std::complex<float>* X,const int incX)
{
	cblas_cscal(N,&alpha,X,incX);
}
void copt_blas_scal(const int N,const std::complex<double>& alpha,std::complex<double>* X,const int incX)
{
	cblas_zscal(N,&alpha,X,incX);
}
void copt_blas_scal(const int N,const float alpha,std::complex<float>* X,const int incX)
{
	cblas_csscal(N,alpha,X,incX);
}
void copt_blas_scal(const int N,const double alpha,std::complex<double>* X,const int incX)
{
	cblas_zdscal(N,alpha,X,incX);
}

/*				summation
*/

//			level 2
/*				multipliation
 */
// void copt_blas_gemm()


}// End of namespace blas

#else
/*					simple blas operations
 */
namespace blas{
template<class eT>
void copt_blas_copy( const int N,const eT* X,const int incX, eT* Y,const int incY)
{
	for ( int i = 0 ; i < N ; ++ i ){
		Y[i*incY] = X[i*incX];
	}
}
template<class eT>
void copt_blas_swap (const int N,eT* X,const int incX,eT* Y,const int incY)
{
	for ( int i = 0 ; i < N ; ++ i ){
		eT t = X[i*incX];
		X[i*incX] = Y[i*incY];
		Y[i*incY] = t;
	}
}
template<class eT>
eT copt_blas_dot(const int N,const eT* X,const int incX,const eT* Y,const int incY)
{
	eT t(0.0);
	for ( int i = 0 ; i < N ; ++ i ){
		t += X[i*incX]*Y[i*incY];
	}
	return t;
}
template<class eT>
void copt_blas_scal(const int N,const eT alpha,eT* X,const int incX)
{
	for (int i = 0 ; i < N ; ++ i ){
		X[i*incX] = X[i*incX]*alpha;
	}
}
}

// End of ifdef CBLAS
#endif


#endif