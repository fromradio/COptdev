//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef UMFPACK_WRAPPER_HPP
#define UMFPACK_WRAPPER_HPP

#ifdef USE_UMFPACK
namespace COPT
{

/** wrappers for umf functions */
//%{

inline void umfpack_free_numeric(void **Numeric,double)
{ 
	umfpack_dl_free_numeric(Numeric); 
	*Numeric = NULL;
}

// inline void umfpack_free_numeric(void **Numeric,std::complex<double>)
// {
// 	umfpack_zl_free_numeric(Numeric);
// 	*Numeric = NULL;
// }

inline void umfpack_free_symbolic(void **Symbolic,double)
{
	umfpack_dl_free_symbolic(Symbolic);
	*Symbolic = NULL;
}

// inline void umfpack_free_symbolic(void **Symbolic,std::complex<double>)
// {
// 	umfpack_zl_free_symbolic(Symbolic);
// 	*Symbolic = NULL;
// }

inline COPTlong umfpack_symbolic(
	COPTlong rows,COPTlong cols,
	const COPTlong* colptr,const COPTlong* rowind,const double* vals,void **Symbolic,
	const double* Control,double* Info)
{
	return umfpack_dl_symbolic(rows,cols,colptr,rowind,vals,Symbolic,Control,Info);
}

inline COPTlong umfpack_symbolic(
	longsize rows,longsize cols,
	const longsize* colptr,const longsize* rowind,const double*vals,void **Symbolic,
	const double* Control , double *Info)
{
	return umfpack_dl_symbolic((COPTlong)rows,(COPTlong)cols,(const COPTlong*)colptr,(const COPTlong*)rowind,vals,Symbolic,Control,Info);
}

// inline int umfpack_symbolic(
// 	int rows,int cols,
// 	const COPTlong* colptr,const COPTlong* rowind,const std::complex<double>* vals,void **Symbolic,
// 	const double* Control,double *Info)
// {
// 	return umfpack_zl_symbolic(rows,cols,colptr,rowind,vals,0,Symbolic,Control,Info);
// }

inline COPTlong umfpack_numeric(
	const COPTlong *colptr, const COPTlong* rowind, const double* vals,
	void *Symbolic, void **Numeric,
	const double* Control,double *Info)
{
	return umfpack_dl_numeric(colptr,rowind,vals,Symbolic,Numeric,Control,Info);
}

inline COPTlong umfpack_numeric(
	const longsize* colptr , const longsize* rowind , const double* vals,
	void *Symbolic, void **Numeric,
	const double* Control,double *Info)
{
	return umfpack_dl_numeric((const COPTlong*)colptr,(const COPTlong*)rowind,vals,Symbolic,Numeric,Control,Info);
}

// inline int umfpack_numeric(
// 	const COPTlong* colptr , const COPTlong* rowind ,const std::complex<double>* vals,
// 	void *Symbolic, void **Numeric,
// 	const double* Control,double *Info)
// {
// 	return umfpack_zl_numeric(colptr,rowind,vals,Symbolic,Numeric,Control,Info);
// }

inline COPTlong umfpack_solve(
	int sys,const COPTlong *colptr,const COPTlong* rowind, const double* vals,
	double* X, const double* B, void *Numeric,
	const double *Control,double *Info)
{
	return umfpack_dl_solve(sys,colptr,rowind,vals,X,B,Numeric,Control,Info);
}

inline COPTlong umfpack_solve(
	int sys,const longsize* colptr, const longsize* rowind, const double* vals,
	double* X, const double* B, void *Numeric,
	const double* Control, double* Info)
{
	return umfpack_dl_solve(sys,(const COPTlong*)colptr,(const COPTlong*)rowind,vals,X,B,Numeric,Control,Info);
}

// inline int umfpack_solve(
// 	int sys,const COPTlong *colptr,const COPTlong* rowind, const double*vals,
// 	std::complex<double>* X,const std::complex<double>* B,void *Numeric,
// 	const double*Control,double *Info)
// {
// 	return umfpack_zl_solve(sys,colptr,rowind,vals,0,X,0,B,0,Numeric,Control,Info);
// }

//%}

/*			A wrapper for UMFPack solving sparse linear system.
 *			The input of the method is a COPT sparse matrix type.
 *			The wrapper first analyze the input matrix and factorize it.
 *			Warning and error information is provided
 */
template<class SpMatrix>
class UMFLinearSolver
	:
	noncopyable
{
private:
	typedef		typename SpMatrix::ScalarType		Scalar;
	typedef 	typename SpMatrix::Size 			Size;

	/**	private variables */
	//%{

	/** the sparse matrix */
	const SpMatrix& 	__mat;

	/** umfpack numeric */
	void*				__numeric;

	/** umfpack numeric */
	void*				__symbolic;

	/** UMFPack Control */
	double* 			__control;

	/** UMFPack Info */
	double*				__info;

	/** whether warning happens */
	bool				__haswarning;

	/** whether error happens */
	bool				__haserror;

	/** status information */
	std::string 		__info_str;

	//%}

	UMFLinearSolver();

	void analyzeStatus(const int status);
public:

	UMFLinearSolver( const SpMatrix& mat);
	~UMFLinearSolver();

	/** given the matrix first analyze the symbolic */
	void analyzeSymbolic();

	/** analyze numeric */
	void analyzeNumeric();

	/** solve a linear system */
	VectorBase<Scalar,Size> solve(const VectorBase<Scalar,Size>& vec);

	/** print information */
	void printInfo();

};

/** 		Implementation			*/

template<class SpMatrix>
void UMFLinearSolver<SpMatrix>::analyzeStatus( const int status )
{
	if ( status == 0 )
	{
		// everythign is ok
		return;
	}
	else if ( status > 0 )
	{
		// there is some warning
		__haswarning = true;
		switch (status)
		{
			case 1:
			{
				__info_str.append("UMFPack Warning: matrix is singular. There are exact zeros on the diagonal of U.\n");
			}
			break;
			case 2:
			{
				__info_str.append("UMFPack Warning: the determinant is nonzero, but smaller in magnitude than the smallest positive floating-point number.\n");
			}
			break;
			case 3:
			{
				__info_str.append("UMFPack Warning: the determinant is larger in magnitude than the largest positive floating-point number (IEEE Inf).\n");
			}
			break;
			default:
			{
				__info_str.append("UMFPack Warning: unknown warning happes!\n");
			}
			break;
		}
	}
	else
	{
		// error happens
		__haserror = true;
		switch (status)
		{
			case -1:
			{
				__info_str.append("UMFPack Error: out of memory!\n");
			}
			break;
			case -3:
			{
				__info_str.append("UMFPack Error: Numeric object is invalid. You should check umfpack user guide for further infomation!\n");
			}
			break;
			case -4:
			{
				__info_str.append("UMFPack Error: Symbolic object is invalid. You should check umfpack user guide for further infomation!\n");
			}
			break;
			case -5:
			{
				__info_str.append("UMFPack Error: A NULL pointer is passed while it need to be present. Check that whether Symbolic or Numeric is NULL or not!\n");
			}
			break;
			case -6:
			{
				__info_str.append("UMFPack Error: the number of rows and columns should be greater than zero!\n");
			}
			break;
			case -8:
			{
				__info_str.append("UMFPack Error: the matrix is invalid!\n");
			}
			break;
			case -11:
			{
				__info_str.append("UMFPack Error: different pattern now. Check that whether that you change the matrix between symbolic and numeric factorization.\n");
			}
			break;
			case -13:
			{
				__info_str.append("UMFPack Error: sys system argument is invalid!\n");
			}
			break;
			case -15:
			{
				__info_str.append("UMFPack Error: the provided permutation vector is invalid!\n");
			}
			break;
			case -17:
			{
				__info_str.append("UMFPack Error: file I/O error happens when saving or loading Numeric or Symbolic object!\n");
			}
			break;
			case -18:
			{
				__info_str.append("UMFPack Error: the ordering method failed!\n");
			}
			break;
			case -911:
			{
				__info_str.append("UMFPack Error: an internal error has occured!\n");
			}
			break;
			default:
			{
				__info_str.append("UMFPack Error: unknown error occured!\n");
			}
			break;
		}
	}
}

template<class SpMatrix>
UMFLinearSolver<SpMatrix>::UMFLinearSolver(
	const SpMatrix& mat)
	:
	__mat(mat),
	__control(new double[UMFPACK_CONTROL]),
	__info(new double[UMFPACK_INFO]),
	__haswarning(false),
	__haserror(false)
{
	analyzeSymbolic();
	analyzeNumeric();

	if(__haswarning||__haserror)
		printInfo();
	else
		__info_str.append("factorization success!");
}

template<class SpMatrix>
UMFLinearSolver<SpMatrix>::~UMFLinearSolver()
{
	umfpack_free_symbolic(&__symbolic,Scalar());
	umfpack_free_numeric(&__numeric,Scalar());
	delete[] __control;
	delete[] __info;
}

template<class SpMatrix>
void UMFLinearSolver<SpMatrix>::analyzeSymbolic()
{
	int status = umfpack_symbolic(__mat.rows(),__mat.cols(),__mat.columnPointer(),__mat.rowIndex(),__mat.values(),&__symbolic,__control,__info);
	analyzeStatus(status);
}

template<class SpMatrix>
void UMFLinearSolver<SpMatrix>::analyzeNumeric()
{
	int status = umfpack_numeric(__mat.columnPointer(),__mat.rowIndex(),__mat.values(),__symbolic,&__numeric,__control,__info);
	analyzeStatus(status);
}

template<class SpMatrix>
VectorBase<typename SpMatrix::ScalarType,typename SpMatrix::Size> UMFLinearSolver<SpMatrix>::solve(const VectorBase<Scalar,Size>& vec)
{
	Scalar* x= new Scalar[vec.size()];
	umfpack_solve(UMFPACK_A,__mat.columnPointer(),__mat.rowIndex(),__mat.values(),x,vec.dataPtr(),__numeric,__control,__info);
	return VectorBase<Scalar,Size>(vec.size(),x);
	delete[]x;
}

template<class SpMatrix>
void UMFLinearSolver<SpMatrix>::printInfo()
{
	std::cerr<<__info_str<<std::endl;
}

}// End of namespace COPT
#endif

#endif