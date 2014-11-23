//		Copyright (C) Songtao Guo, guost@mathu.cn
//		Copyright (C) MathU
//			Code has been reviewed and modified by Ruimin Wang, ruimin.wang13@gmail.com, wangrm@mathu.cn

#ifndef LEAST_SQUARE_H
#define LEAST_SQUARE_H


namespace COPT
{
/*
*            Currently,three methods are contained:
*            Least Mean Square Method
*            Least Square Method
*            Recursive Least Square Method
*/
template<class Scalar>
class LeastSquaresSolver
{
public:
	enum SolverType{
		LMS,       //Least Mean Square method
		LS,        //Least Square method
		RLS,       //Recursive Least Square method
	};
	LeastSquaresSolver(
		const MatrixBase<Scalar>& A,
		const VectorBase<Scalar>& b,
		const Scalar mu = 0.01,
		const Scalar lam = 1,
		const Scalar delta = 250,
		const SolverType type = LMS
		);

	// set the type of the solver:LMS,LS or RLS
	void setType(const SolverType type);

    // interface for solving the problem
	void solve(VectorBase<Scalar>& x);

	// print the solver information
	void printInfo();

private:
	// the coefficient matrix
	const MatrixBase<Scalar>&   __A;
	// constant vector
	const VectorBase<Scalar>&   __b;
	// parameter of Least Mean Square method
	Scalar                      __mu;
	// parameter of Recursive Least Square method
	Scalar                      __lam;
	// parameter of Recursive Least Square method
	Scalar                      __delta;
	// the final result
	VectorBase<Scalar>          __x;
	//the type of the solver,LMS as default
	SolverType                  __type;

    //the main algorithm
	void            doSolve();

public:
	void LeastMeanSquareMethod(
		const MatrixBase<Scalar>& A,
		const VectorBase<Scalar>& b,
		const Scalar mu,
		VectorBase<Scalar>& x
		);
	void LeastSquareMethod(
		const MatrixBase<Scalar>& A,
		const VectorBase<Scalar>& b,
		VectorBase<Scalar>& x
		);
	void RLS_Method(
		const MatrixBase<Scalar>& A,
		const VectorBase<Scalar>& b,
		VectorBase<Scalar>& x,
		const Scalar lam = 1,
		const Scalar delta = 250);

private:
	void LeastMeanSquareUpdate(
	const VectorBase<Scalar>& a,
	const Scalar b,
	const Scalar mu,
	VectorBase<Scalar>& x
	);
	void RLS_learning(
	const VectorBase<Scalar>& a,
	const Scalar b,
	VectorBase<Scalar>& x,
	MatrixBase<Scalar>& p,
	const Scalar lam,
	const Scalar delta
	);
};


template<class Scalar>
LeastSquaresSolver<Scalar>::LeastSquaresSolver(
	const MatrixBase<Scalar>& A,
	const VectorBase<Scalar>& b,
	const Scalar mu,
	const Scalar lam,
	const Scalar delta,
	const SolverType type
	)
    :
    __A(A),
    __b(b),
    __mu(mu),
    __lam(lam),
    __delta(delta),
    __type(type)	
{
}

/*
*        set the type of the solver: LMS , LS or RLS
*/
template<class Scalar>
void LeastSquaresSolver<Scalar>::setType(const SolverType type)
{
	__type = type;
}


/*      Kernel function
*      Solving the problem
*/
template<class Scalar>
void LeastSquaresSolver<Scalar>::doSolve()
{
	switch(__type)
	{
	case LMS:
	{
		LeastMeanSquareMethod(
			__A,
			__b,
			__mu,
			__x);
	}
	break;
	case LS:
	{
		LeastSquareMethod(
			__A,
			__b,
			__x);
	}
	break;
	case RLS:
	{
		RLS_Method(
			__A,
			__b,
			__x,
			__lam,
			__delta);
	}
	break;
	}
}


/*        Solver interface for users
*         /param x:              inital point
*/
template<class Scalar>
void LeastSquaresSolver<Scalar>::solve(
	VectorBase<Scalar>& x
	)
{
	__x = x;
	doSolve();
}

template<class Scalar>
void LeastSquaresSolver<Scalar>::printInfo()
{
	std::cout<<"The result of LeastSquares solver by ";
	switch(__type){
	case LMS:
	{
		std::cout<<"Least Mean Square Method";
	}
	break;
	case LS:
	{
		std::cout<<"Least Square Method";
	}
	break;
	case RLS:
	{
		std::cout<<"Recursive Least Square Method";
	}
	break;
	}
	std::cout<<":"<<std::endl;
	std::cout<<__x<<std::endl;
}

#endif






/*
*    Least Mean Square algorithm update step
*    /param a:         row of the coefficient matrix
*    /param b:         target value corresponding to a
*    /param mu:        step size
*    /param x:         weight we find after each step
*/
template<class Scalar>
void LeastSquaresSolver<Scalar>::LeastMeanSquareUpdate(
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
void LeastSquaresSolver<Scalar>::LeastMeanSquareMethod(
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
void LeastSquaresSolver<Scalar>::LeastSquareMethod(
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
void LeastSquaresSolver<Scalar>::RLS_learning(
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
void LeastSquaresSolver<Scalar>::RLS_Method(
	const MatrixBase<Scalar>& A,
	const VectorBase<Scalar>& b,
	VectorBase<Scalar>& x,
	const Scalar lam,
	const Scalar delta
	)
{
	int n = A.cols();
	MatrixBase<Scalar> p = delta*MatrixBase<Scalar>::identity(n,n);
	for(int i = 0;i < n;i++)
		RLS_learning(A.row(i),b[i],x,p,lam,delta);
}

}
