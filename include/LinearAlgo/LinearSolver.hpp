//		Copyright (C) ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef LINEAR_SOLVER_HPP__
#define LINEAR_SOLVER_HPP__


namespace COPT
{
/*			A general design for solving linear system 'LinearSolver'.
 *			LinearSolver takes the type of matrix as template. The class
 *			stores the size and lda of input matrix.
 *			A derived class from LinearSolver has to implement several
 *			special functions as the computation of the linear system and
 *			the solving of the system.
 */
template<class Matrix>
class LinearSolver
	:
	COPTObject,
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

	/** virtual function for implementation */
	//%{
	/** computation of the solver */
	virtual void doCompute( const Matrix& mat ) = 0;
	/** solve a linear system with a right hand vectors */
	virtual Vector doSolve( const Vector& b ) = 0;
	virtual Matrix doSolve( const Matrix& b ) = 0;
	//%}

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
	/** lda of the matrix */
	void setLDA( const index lda );
	index lda( ) const;
	/** row number */
	void setRowNum( const index m );
	index rowNum() const;
	/** column number */
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

/* 
 *		some construction of COPT linear solvers 
 */
template<class Matrix> class LU;
template<class Matrix> class QR;
template<class Matrix> class Cholesky;
template<class Matrix> class EigenSolver;

void setLinearSolver()
{

}

}
#endif