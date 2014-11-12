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
	const Matrix& A,
	const Vector& b,
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

//Recursive Least Square
/*
*    Recursive Least Square algorithm update step
*    /param a:         row of the coefficient matrix
*    /param b:         target value corresponding to a
*    /param x:         weight we find after each step   
*    /param p:         a trival variable    
*    /param lam:       related to the forgetting factor
*    /param delta:     in order to define P(0)   
*/
template<class Vector,class Matrix>
void RLS_learning(
	const Vector& a,
	const typename Vector::ScalarType b,
	Vector& x,
	Matrix& p,
	const typename Vector::ScalarType lam,
	const typename Vector::ScalarType delta
	)
{
	Vector pai = p.transpose()*a;
	typename Vector::ScalarType gama = lam + pai.dot(a);
	Vector k = pai*(1.0/gama);
	typename Vector::ScalarType alpha = b - x.dot(a);
	x = x + k*alpha;
	Matrix pp = k.mulTrans(pai); 
	p = (1.0/lam)*(p - pp);
}
/*
*    Recursive Least Square algorithm 
*    /param A:         the coefficient matrix
*    /param b:         constant vector
*    /param x:         weight we want to find   
*    /param lam:       related to the forgetting factor
*    /param delta:     in order to define P(0)   
*/
template<class Vector,class Matrix>
void RLS_Method(
	const Matrix& A,
	const Vector& b,
	Vector& x,
	const typename Vector::ScalarType lam = 1,
	const typename Vector::ScalarType delta = 250
	)
{
	int n = A.cols();
	Matrix p = delta*Matrix::identity(n,n);
	for(int i = 0;i < n;i++)
		RLS_learning(A.row(i),b[i],x,p,lam,delta);
}

}


#endif