//		Copyright (C) Songtao Guo, guost@mathu.cn
//		Copyright (C) MathU
//			Code has been reviewed and modified by Ruimin Wang, ruimin.wang13@gmail.com, wangrm@mathu.cn

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
template<class Scalar>
void LeastMeanSquareUpdate(
	const VectorBase<Scalar>& a,
	const Scalar b,
	const Scalar mu,
	VectorBase<Scalar>& x
	)
{
	Scalar e = b - a.dot(x);
	x = x + mu*e*a;
}

/*
*    Least Mean Square algorithm
*    /param A:          the coefficient matrix
*    /param b:          constant vector
*    /param mu:         step size
*    /param x:          weight we want to find
*/
template<class Scalar>
void LeastMeanSquareMethod(
	const MatrixBase<Scalar>& A,
	const VectorBase<Scalar>& b,
	const Scalar mu,
	VectorBase<Scalar>& x 
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
template<class Scalar>
void LeastSquareMethod(
	const MatrixBase<Scalar>& A,
	const VectorBase<Scalar>& b,
	VectorBase<Scalar>& x
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
template<class Scalar>
void RLS_learning(
	const VectorBase<Scalar>& a,
	const Scalar b,
	VectorBase<Scalar>& x,
	MatrixBase<Scalar>& p,
	const Scalar lam,
	const Scalar delta
	)
{
	VectorBase<Scalar> pai = p.transpose()*a;
	Scalar gama = lam + pai.dot(a);
	VectorBase<Scalar> k = pai*(1.0/gama);
	Scalar alpha = b - x.dot(a);
	x = x + k*alpha;
	MatrixBase<Scalar> pp = k.mulTrans(pai); 
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
template<class Scalar>
void RLS_Method(
	const MatrixBase<Scalar>& A,
	const VectorBase<Scalar>& b,
	VectorBase<Scalar>& x,
	const Scalar lam = 1,
	const Scalar delta = 250
	)
{
	int n = A.cols();
	MatrixBase<Scalar> p = delta*MatrixBase<Scalar>::identity(n,n);
	for(int i = 0;i < n;i++)
		RLS_learning(A.row(i),b[i],x,p,lam,delta);
}

}


#endif