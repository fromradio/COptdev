// This file is part of COPT, a light-weight C++ based optimization open source library
//
// Copyright (C) 2015 Songtao Guo <guost@mathu.cn>
// Copyright (C) 2015 MathU
//			Reviewed by Ruimin Wang <ruimin.wang13@gmail.com>, <wangrm@mathu.cn>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.



#ifndef LEAST_SQUARE_HPP__
#define LEAST_SQUARE_HPP__


namespace COPT
{
/*         Least Mean Square solver for LeastSquares problem
*
*
*/ 
template<class Problem,class Time = NoTimeStatistics>
class LeastMeanSquareSolver
    :
    public GeneralSolver<typename Problem::KernelTrait,Time>
{
private:
	typedef typename Problem::KernelTrait::index		index;
	typedef typename Problem::KernelTrait::scalar 		scalar;
	typedef typename Problem::KernelTrait::Vector 		Vector;
	typedef typename Problem::KernelTrait::Matrix 		Matrix;
	/** private variables */
	//%{
	/** the reference to the problem */
	const Problem&              __p;

	/** the coefficient matrix */
	Matrix                      __A;

	/** constant vector */
	Vector                      __b;

	/** iteration number */
	scalar                      __iter_n;

	/** rows of __A */
	Vector                      __Ar;

	/** elements of __b */
	scalar                      __be;

	/** parameter of Least Mean Square method */
	scalar                      __mu;

	/** the final result */
	Vector                      __x;

	//%} end of variables

	// void init();

	LeastMeanSquareSolver();

	void doSolve();

	void solvingBegin();
	void doCompute();
	scalar doOneIteration();


public:

//	typedef leastmeansquare_solver              ObjectCategory;
	/** constructor and deconstructor */
	//%{
	LeastMeanSquareSolver(
		const Problem& s ,
		const scalar mu = 0.01);

	//%} end of constructor and deconstructor

	/** setter and getter */
	//%{
	/** the final result of proximal solver */
	const Vector& result() const;
	//%}

	scalar objective() const;

};




/*         Least Square solver for LeastSquares problem
*
*
*/ 

template<class Problem,class Time = NoTimeStatistics>
class LeastSquareSolver
    :
    public GeneralSolver<typename Problem::KernelTrait,Time>
{
private:
	typedef typename Problem::KernelTrait::index		index;
	typedef typename Problem::KernelTrait::scalar 		scalar;
	typedef typename Problem::KernelTrait::Vector 		Vector;
	typedef typename Problem::KernelTrait::Matrix 		Matrix;
	/** private variables */
	//%{
	/** the reference to the problem */
	const Problem&              __p;

	/** the coefficient matrix */
	Matrix                      __A;

	/** constant vector */
	Vector                      __b;

	/** the final result */
	Vector                      __x;

	//%} end of variables

	// void init();

	LeastSquareSolver();

	void doSolve();

	void solvingBegin();
	void doCompute();
	scalar doOneIteration();
    scalar objective() const;


public:

//	typedef leastsquare_solver              ObjectCategory;
	/** constructor and deconstructor */
	//%{
	LeastSquareSolver(
		const Problem& s);

	//%} end of constructor and deconstructor

	/** setter and getter */
	//%{
	/** the final result of proximal solver */
	const Vector& result() const;
	//%}

};




/*         Recursive Least Square solver for LeastSquares problem
*
*
*/ 

template<class Problem,class Time = NoTimeStatistics>
class RecursiveLeastSquareSolver
    :
    public GeneralSolver<typename Problem::KernelTrait,Time>
{
private:
	typedef typename Problem::KernelTrait::index		index;
	typedef typename Problem::KernelTrait::scalar 		scalar;
	typedef typename Problem::KernelTrait::Vector 		Vector;
	typedef typename Problem::KernelTrait::Matrix 		Matrix;
	/** private variables */
	//%{
	/** the reference to the problem */
	const Problem&              __p;

	/** the coefficient matrix */
	Matrix                      __A;

	/** constant vector */
	Vector                      __b;

	/** iteration number , also rows of A */
	scalar                      __iter_n;

	/** cols number of A */
	scalar                      __iter_m;

	/** rows of __A */
	Vector                      __Ar;

	/** elements of __b */
	scalar                      __be;

	/** parameter of Recursive Least Square method */
	scalar                      __lam;

	/** parameter of Recursive Least Square method */
	scalar                      __delta;

	/** intermediate variable */
	Matrix                      __inter_p;

	/** the final result */
	Vector                      __x;

	//%} end of variables

	// void init();

	RecursiveLeastSquareSolver();

	void doSolve();

	void solvingBegin();
	void doCompute();
	scalar doOneIteration();
	scalar objective() const;


public:

//	typedef recursiveleastsquare_solver              ObjectCategory;
	/** constructor and deconstructor */
	//%{
	RecursiveLeastSquareSolver(
		const Problem& s ,
		const scalar lam = 1 ,
		const scalar delta = 250);

	//%} end of constructor and deconstructor

	/** setter and getter */
	//%{
	/** the final result of proximal solver */
	const Vector& result() const;
	//%}

};



	
/*		The Least Squares problem class
 *
 */
template<class kernel>
class LeastSquaresProblem
	:
	public VectorProblem<kernel>
{
private:
	typedef typename kernel::scalar 		scalar;
	typedef typename kernel::Matrix 		Matrix;
	typedef typename kernel::Vector 		Vector;

	const Matrix& 				__A;
	const Vector& 				__b;

public:

	typedef kernel							KernelTrait;
	// typedef leastsquares_problem					ObjectCategory;

	LeastSquaresProblem(
		const Matrix& A,
		const Vector& b);

	// void leastMeanSquareSolve();
	// void leastSquareSolve();
	// void recursiveLeastSquareSolve();


	/** check whether the input is valid */
	bool isValidInput( const Vector& x ) const;

	/** check whether the problem is a valid problem */
	bool isValid( ) const;


	/** getter and setter */
	//%{
	const Matrix& matA() const;
	const Vector& obB() const;
	scalar objective( const Vector& x) const;
	//%}
};


/*********************Implementation of LeastMeanSquareSolver ******************/


template<class Problem,class Time>
void LeastMeanSquareSolver<Problem,Time>::doSolve()
{
	this->__iter_num = 0;
	do
	{
		this->__estimated_error = this->oneIteration();
		if(++this->__iter_num>=__iter_n)
		{
			this->__terminal_type = this->MaxIteration;
			break;
		}
	}while(this->__estimated_error>this->__thresh);
	if(this->__estimated_error<=this->__thresh)
		this->__terminal_type = this->Optimal;
}

template<class Problem,class Time>
void LeastMeanSquareSolver<Problem,Time>::doCompute()
{
	__A = __p.matA();
	__b = __p.obB();
}

template<class Problem,class Time>
LeastMeanSquareSolver<Problem,Time>::LeastMeanSquareSolver( 
	const Problem& s , 
	const scalar mu)
	:
	__p(s),
	__mu(mu)
{
	this->doCompute();
}

template<class Problem,class Time>
void LeastMeanSquareSolver<Problem,Time>::solvingBegin()
{
	__x.resize(__A.cols());
	__iter_n = __A.rows();
}

template<class Problem,class Time>
typename LeastMeanSquareSolver<Problem,Time>::scalar LeastMeanSquareSolver<Problem,Time>::objective() const
{
	return __p.objective(__x);
}

template<class Problem,class Time>
typename LeastMeanSquareSolver<Problem,Time>::scalar LeastMeanSquareSolver<Problem,Time>::doOneIteration()
{
	Vector xp = __x;
	__Ar = __A.row(this->__iter_num);
	__be = __b[this->__iter_num];
	scalar e = __be - __Ar.dot(__x);
	__x = __x + __mu*e*__Ar;
	return (__A*__x-__b).squaredNorm();
}

template<class Problem,class Time>
const typename LeastMeanSquareSolver<Problem,Time>::Vector& LeastMeanSquareSolver<Problem,Time>::result() const
{
	return __x;
}



/*********************Implementation of LeastSquareSolver ******************/

template<class Problem,class Time>
void LeastSquareSolver<Problem,Time>::doSolve()
{
	this->__estimated_error = this->oneIteration();
	if(this->__estimated_error<=this->__thresh)
		this->__terminal_type = this->Optimal;
}

template<class Problem,class Time>
void LeastSquareSolver<Problem,Time>::doCompute()
{
	__A = __p.matA();
	__b = __p.obB();
}

template<class Problem,class Time>
LeastSquareSolver<Problem,Time>::LeastSquareSolver( 
	const Problem& s)
	:
	__p(s)
{
	this->doCompute();
}

template<class Problem,class Time>
void LeastSquareSolver<Problem,Time>::solvingBegin() 
{
	__x.resize(__A.cols());
}

template<class Problem,class Time>
typename LeastSquareSolver<Problem,Time>::scalar LeastSquareSolver<Problem,Time>::objective() const
{
	return __p.objective(__x);
}

template<class Problem,class Time>
typename LeastSquareSolver<Problem,Time>::scalar LeastSquareSolver<Problem,Time>::doOneIteration()
{
	Vector xp = __x;
	__x = (__A.transpose()*__A).solve(__A.transpose()*__b);
	return (__A*__x-__b).squaredNorm();
}

template<class Problem,class Time>
const typename LeastSquareSolver<Problem,Time>::Vector& LeastSquareSolver<Problem,Time>::result() const
{
	return __x;
}



/*********************Implementation of RecursiveLeastSquareSolver ******************/

template<class Problem,class Time>
void RecursiveLeastSquareSolver<Problem,Time>::doSolve()
{
	this->__iter_num = 0;
	do
	{
		this->__estimated_error = this->oneIteration();
		if(++this->__iter_num >= __iter_n)
		{
			this->__terminal_type = this->MaxIteration;
			break;
		}
	}while(this->__estimated_error>this->__thresh);
	if(this->__estimated_error<=this->__thresh)
		this->__terminal_type = this->Optimal;
}

template<class Problem,class Time>
void RecursiveLeastSquareSolver<Problem,Time>::doCompute()
{
	__A = __p.matA();
	__b = __p.obB();
}

template<class Problem,class Time>
RecursiveLeastSquareSolver<Problem,Time>::RecursiveLeastSquareSolver( 
	const Problem& s , 
	const scalar lam ,
	const scalar delta)
	:
	__p(s),
	__lam(lam),
	__delta(delta)
{
	this->doCompute();
}

template<class Problem,class Time>
void RecursiveLeastSquareSolver<Problem,Time>::solvingBegin() 
{
	__x.resize(__A.cols());
	__iter_n = __A.rows();
	__iter_m = __A.cols();
	__inter_p = __delta*Matrix::identity(__iter_m,__iter_m);
}

template<class Problem,class Time>
typename RecursiveLeastSquareSolver<Problem,Time>::scalar RecursiveLeastSquareSolver<Problem,Time>::objective() const
{
	return __p.objective(__x);
}

template<class Problem,class Time>
typename RecursiveLeastSquareSolver<Problem,Time>::scalar RecursiveLeastSquareSolver<Problem,Time>::doOneIteration()
{
	Vector xp = __x;
	__Ar = __A.row(this->__iter_num);
	__be = __b[this->__iter_num];
	Vector pai = __inter_p.transpose()*__Ar;
	scalar gama = __lam + pai.dot(__Ar);
	Vector k = pai*(1.0/gama);
	scalar alpha = __be - __x.dot(__Ar);
	__x = __x + k*alpha;
	Matrix pp = k.mulTrans(pai); 
	__inter_p = (1.0/__lam)*(__inter_p - pp);
	return (__A*__x-__b).squaredNorm();
}

template<class Problem,class Time>
const typename RecursiveLeastSquareSolver<Problem,Time>::Vector& RecursiveLeastSquareSolver<Problem,Time>::result() const
{
	return __x;
}



/***********************Implementation of LeastSquaresProblem************************/

template<class kernel>
LeastSquaresProblem<kernel>::LeastSquaresProblem(
	const Matrix& A,
	const Vector& b)
	:
	__A(A),
	__b(b)
{  
}

// template<class kernel>
// void LeastSquaresProblem<kernel>::leastMeanSquareSolve()
// {
// }

// template<class kernel>
// void LeastSquaresProblem<kernel>::leastSquareSolve()
// {
// }

// template<class kernel>
// void LeastSquaresProblem<kernel>::recursiveLeastSquareSolve()
// {
// }

template<class kernel>
bool LeastSquaresProblem<kernel>::isValid( ) const
{
	if (__A.rows()!=__b.size())
		return false;
	else 
		return true;
}

template<class kernel>
bool LeastSquaresProblem<kernel>::isValidInput( const Vector& x ) const
{
	if (__A.cols() != x.size() )
		return false;
	else 
		return true;
}

template<class kernel>
const typename kernel::Matrix& LeastSquaresProblem<kernel>::matA() const
{
	return __A;
}

template<class kernel>
const typename kernel::Vector& LeastSquaresProblem<kernel>::obB() const
{
	return __b;
}

template<class kernel>
typename kernel::scalar LeastSquaresProblem<kernel>::objective( const Vector& x ) const
{
	return (__A*x-__b).squaredNorm();
}
}


#endif
