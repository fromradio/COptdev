//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef EIGEN_WRAPPER_HPP__
#define EIGEN_WRAPPER_HPP__


namespace COPT
{

/*			A wrapper from open source library 'Eigen's linear solver.
 *			The ldlt solver is chosen since it is the fastest one due to
 *			the description of 'Eigen'.
 */
template<class Matrix>
class EigenSolver
	:
	public LinearSolver<Matrix>
{
private:
	typedef typename Matrix::scalar 				scalar;
	typedef typename Matrix::index 					index;
	typedef VectorBase<scalar,index>				Vector;


	void doCompute(const Matrix& mat){}
	Vector doSolve( const Vector& b ){return Vector();}
	Matrix doSolve( const Matrix& b ){return Matrix();}

public:
	
	EigenSolver(){}
	~EigenSolver(){}

	void clear(){}
};


}

#endif