//		Copyright (C) Songtao Guo, guost@mathu.cn
//		Copyright (C) MathU

#ifndef LEAST_SQUARE_H
#define LEAST_SQUARE_H

// Least mean square

namespace COPT
{

/*
*    Least Mean Square algorithm update step
*    /param a:         row of the coefficient matrix
*    /param b:         target value corresponding to a
*    /param mu:        step size
*    /param x:         weight we find after each step
*/
template<class Vector>
void LeastMeanSquareUpdate(
	const Vector& a,
	const typename Vector::ScalarType b,
	const typename Vector::ScalarType mu,
	Vector& x
	)
{
	typename Vector::ScalarType e = b - a.dot(x);
	x = x + mu*e*a;
}

/*
*    Least Mean Square algorithm
*    /param A:          the coefficient matrix
*    /param b:          constant vector
*    /param mu:         step size
*    /param x:          weight we want to find
*/
template<class Vector,class Matrix>
void LeastMeanSquareMethod(
	const Matrix A,
	const Vector b,
	const typename Vector::ScalarType mu,
	Vector& x 
	)
{
	int n = A.rows();
	for(int i = 0;i < n;i++)
		LeastMeanSquareUpdate(A.row(i),b[i],mu,x);

}

/*
*    Least Square Method
*    /param A:          the coefficient matrix
*    /param b:          constant vector
*    /param x:          weight we want to find
*/
template<class Vector,class Matrix>
void LeastSquareMethod(
	const Matrix A,
	const Vector b,
	Vector& x
	)
{
	x = (A.transpose()*A).solve(A.transpose()*b);
}

}


#endif