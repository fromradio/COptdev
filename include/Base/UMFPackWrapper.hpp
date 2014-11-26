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

inline int umfpack_symbolic(
	int rows,int cols,
	const COPTlong* colptr,const COPTlong* rowind,const double* vals,void **Symbolic,
	const double* Control,double* Info)
{
	return umfpack_dl_symbolic(rows,cols,colptr,rowind,vals,Symbolic,Control,Info);
}

// inline int umfpack_symbolic(
// 	int rows,int cols,
// 	const COPTlong* colptr,const COPTlong* rowind,const std::complex<double>* vals,void **Symbolic,
// 	const double* Control,double *Info)
// {
// 	return umfpack_zl_symbolic(rows,cols,colptr,rowind,vals,0,Symbolic,Control,Info);
// }

inline int umfpack_numeric(
	const COPTlong *colptr, const COPTlong* rowind, const double* vals,
	void *Symbolic, void **Numeric,
	const double* Control,double *Info)
{
	return umfpack_dl_numeric(colptr,rowind,vals,Symbolic,Numeric,Control,Info);
}

// inline int umfpack_numeric(
// 	const COPTlong* colptr , const COPTlong* rowind ,const std::complex<double>* vals,
// 	void *Symbolic, void **Numeric,
// 	const double* Control,double *Info)
// {
// 	return umfpack_zl_numeric(colptr,rowind,vals,Symbolic,Numeric,Control,Info);
// }

inline int umfpack_solve(
	int sys,const COPTlong *colptr,const COPTlong* rowind, const double* vals,
	double* X, const double* B, void *Numeric,
	const double *Control,double *Info)
{
	return umfpack_dl_solve(sys,colptr,rowind,vals,X,B,Numeric,Control,Info);
}

// inline int umfpack_solve(
// 	int sys,const COPTlong *colptr,const COPTlong* rowind, const double*vals,
// 	std::complex<double>* X,const std::complex<double>* B,void *Numeric,
// 	const double*Control,double *Info)
// {
// 	return umfpack_zl_solve(sys,colptr,rowind,vals,0,X,0,B,0,Numeric,Control,Info);
// }

//%}

template<class Scalar,class Size = size_t >
class UMFLinearSolver
	:
	noncopyable
{
private:

	// typedef 		integer_type<Size>		Integer;

	/**	private variables */
	//%{

	/** number of rows */
	const Size 		__rows;

	/** number of columns */
	const Size 		__cols;

	/** number of non-zero elements */
	const Size  	__nnz;

	/** column pointer */
	const Size* 	__colptr;

	/** row indices */
	const Size* 	__rowind;

	/** pointer to values */
	const Scalar* 	__vals;

	/** umfpack numeric */
	void**			__numeric;

	/** umfpack numeric */
	void**			__symbolic;

	/** UMFPack Control */
	double* 		__control;

	/** UMFPack Info */
	double*			__info;

	//%}

	UMFLinearSolver();
public:

	UMFLinearSolver(
		const Size rows,
		const Size cols,
		const Size nnz,
		const Size* colptr,
		const Size* rowind,
		const Scalar* vals);
	~UMFLinearSolver();

	/** given the matrix first analyze the symbolic */
	
};

/** 		Implementation			*/
template<class Scalar,class Size>
UMFLinearSolver<Scalar,Size>::UMFLinearSolver(
	const Size rows,
	const Size cols,
	const Size nnz,
	const Size* colptr,
	const Size* rowind,
	const Scalar* vals)
	:
	__rows(rows),
	__cols(cols),
	__nnz(nnz),
	__colptr(colptr),
	__rowind(rowind),
	__vals(vals),
	__control(new double[UMFPACK_CONTROL]),
	__info(new double[UMFPACK_INFO])
{
}

template<class Scalar,class Size>
UMFLinearSolver<Scalar,Size>::~UMFLinearSolver()
{
	umfpack_free_symbolic(__symbolic,Scalar(),Size());
	umfpack_free_numeric(__numeric,Scalar(),Size());
	delete[] __control;
	delete[] __info;
}

}// End of namespace COPT
#endif

#endif