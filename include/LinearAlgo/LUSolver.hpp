//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef LU_SOLVER_HPP__
#define LU_SOLVER_HPP__

namespace COPT
{

template<class Matrix>
class LinearSolver
	:
	noncopyable
{

	typedef typename Matrix::scalar 		scalar;
	typedef typename Matrix::index 			index;
	typedef VectorBase<scalar,index>		Vector;

	/** private varibles */
	//%{
	/** the row number of matrix */
	index 					__m;
	/** the column number of matrix */
	index 					__n;
	/** lda of matrix */
	index 					__lda;
	//%}

	virtual void doCompute( const Matrix& mat ) = 0;
	virtual Vector doSolve( const Vector& b ) = 0;
	virtual Matrix doSolve( const Matrix& b ) = 0;

public:

	virtual ~LinearSolver() {}

	/** compute the solver */
	virtual void compute ( const Matrix& mat );

	/** solving one right hand vector */
	virtual Vector solve ( const Vector& b );

	/** solving several right hand vectors */
	virtual Matrix solve ( const Matrix& b );

	/** clear the solver */
	virtual void clear( ) = 0;

	/** square validation */
	virtual void squareValidation() const;

	/** setter and getter */
	//%{
	void setLDA( const index lda );
	index lda( ) const;
	void setRowNum( const index m );
	index rowNum() const;
	void setColNum( const index n );
	index colNum() const;
	//%}
};

/***********Implementation of class 'LinearSolver'****************/

template<class Matrix>
void LinearSolver<Matrix>::compute( const Matrix& mat )
{
	this->doCompute(mat);
}

template<class Matrix>
typename LinearSolver<Matrix>::Vector LinearSolver<Matrix>::solve ( const Vector& b )
{
	return this->doSolve( b );
}

template<class Matrix>
Matrix LinearSolver<Matrix>::solve ( const Matrix& b )
{
	return this->doSolve(b);
}

template<class Matrix>
void LinearSolver<Matrix>::squareValidation() const
{
	if ( __n != __m )
		throw COException("Linear system solving error: the matrix is not square! ");
}

template<class Matrix>
void LinearSolver<Matrix>::setLDA( const index lda )
{
	__lda = lda;
}

template<class Matrix>
typename LinearSolver<Matrix>::index LinearSolver<Matrix>::lda() const
{
	return __lda;
}

template<class Matrix>
void LinearSolver<Matrix>::setRowNum( const index m )
{
	__m = m;
}

template<class Matrix>
typename LinearSolver<Matrix>::index LinearSolver<Matrix>::rowNum() const
{
	return __m;
}

template<class Matrix>
void LinearSolver<Matrix>::setColNum( const index n )
{
	__n = n;
}

template<class Matrix>
typename LinearSolver<Matrix>::index LinearSolver<Matrix>::colNum() const
{
	return __n;
}


/////////////End of implementation of class 'LinearSolver'

template<class Matrix>
class LU
{
private:
	typedef typename Matrix::scalar 		scalar;
	typedef typename Matrix::index 			index;
	typedef VectorBase<scalar,index>		Vector;

	/** the array */
	scalar 						*__a;
	/** the size of array */
	index 						__size;
	/** lda of a */
	index 						__lda;
	/** the row number of the input matrix */
	index 						__rows;
	/** the column number of the input matrix */
	index 						__cols;
	/** piv */
	index 						*__piv;
	/** the factorization information */
	index 						__info;

public:

	/** constructor and deconstructor */
	//%{
	LU ();
	LU ( const Matrix& mat );
	~LU();
	//%}

	/** computation of the linear system */
	void compute( const Matrix& mat );

	/** square validation before solving */
	void squareValidation() const;

	/** solve one right hand vector */
	Vector solve ( const Vector& b );

	/** solve several right hand vectors */
	Matrix solve ( const Matrix& b );

	/** inverse matrix */
	Matrix inverse ( );

	/** clear */
	void clear( );

};

/*************Implementation of class 'LU'*******************/

template<class Matrix>
LU<Matrix>::LU()
	:
	__a(NULL),
	__piv(NULL)
{
}

template<class Matrix>
LU<Matrix>::LU(const Matrix& mat)
	:
	__a(NULL),
	__piv(NULL)
{
	this->compute(mat);
}

template<class Matrix>
LU<Matrix>::~LU()
{
	clear();
}

template<class Matrix>
void LU<Matrix>::compute( const Matrix& mat )
{
	SAFE_DELETE_ARRAY(__a);
	SAFE_DELETE_ARRAY(__piv);
	__a = new scalar[mat.size()];
	__piv = new index[std::min(mat.rows(),mat.cols())];
	__lda = mat.rows();
	__size = mat.size();
	__rows = mat.rows();
	__cols = mat.cols();
	blas::copt_blas_copy(mat.size(),mat.dataPtr(),1,__a,1);
	copt_lapack_getrf(mat.rows(),mat.cols(),__a,mat.rows(),__piv,&__info);
}

template<class Matrix>
void LU<Matrix>::squareValidation() const
{
	if (__rows!=__cols)
		throw COException("Linear system solving error: the matrix is not square! ");
}

template<class Matrix>
typename LU<Matrix>::Vector LU<Matrix>::solve( const Vector& b )
{
	squareValidation();
	if ( __rows != b.size() )
	{
		std::cerr<<"the order of matrix is "<<__rows<<" and the size of vector is "<<b.size()<<std::endl;
		throw COException("Linear system solving error: the size is not consistent!");
	}
	Vector result(b);
	copt_lapack_getrs('N',__rows,1,__a,__lda,__piv,result.dataPtr(),result.size(),&__info);
	return result;
}

template<class Matrix>
Matrix LU<Matrix>::solve( const Matrix& b )
{
	squareValidation();
	if( __rows != b.rows() )
	{
		std::cerr<<"the order of matrix is "<<__rows<<" and the size of right hand vectors are "<<b.rows()<<std::endl;
		throw COException("Linear system solving error: the size is not consistent!");
	}
	Matrix result(b);
	copt_lapack_getrs('N',__rows,b.cols(),__a,__lda,__piv,result.dataPtr(),result.rows(),&__info);
	return result;
}

template<class Matrix>
Matrix LU<Matrix>::inverse( )
{
	squareValidation();
	Matrix result(__rows,__cols,__a);
	copt_lapack_getri(__rows,result.dataPtr(),__rows,__piv,&__info);
	return result;
}

template<class Matrix>
void LU<Matrix>::clear()
{
	SAFE_DELETE_ARRAY(__a);
	SAFE_DELETE_ARRAY(__piv);
}

}

#endif