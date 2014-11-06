#ifndef CBLAS_WRAPPER_H
#define CBLAS_WRAPPER_H

#ifdef COPT_USE_BLAS
#endif
namespace blas{
template<class eT>
void copt_blas_copy( const int N,eT* X,const int incX, eT* Y)
{
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

}// End of namespace blas

#endif