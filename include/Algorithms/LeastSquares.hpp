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

namespace COPT {
/*	"""     'LeastMeanSquareSolver' is one of the solvers of least squares problems. 
 'LeastMeanSquareSolver' solves the least squares problems using least mean squares algorithm which is a stochastic gradient descent method.
 This class derives directly from another class called 'GeneralSolver' and takes 'Problem' and 'Time' as template.
 """ */
template<class Problem, class Time = NoTimeStatistics>
class LeastMeanSquareSolver: public GeneralSolver<typename Problem::KernelTrait,
		Time> {
private:
	/*** the type of index used in 'LeastMeanSquareSolver' ***/
	typedef typename Problem::KernelTrait::index index;
	/*** the type of scalar used in 'LeastMeanSquareSolver' ***/
	typedef typename Problem::KernelTrait::scalar scalar;
	/*** the type of vector used in 'LeastMeanSquareSolver' ***/
	typedef typename Problem::KernelTrait::Vector Vector;
	/*** the type of matrix used in 'LeastMeanSquareSolver' ***/
	typedef typename Problem::KernelTrait::Matrix Matrix;

	/** private variables */
	//%{
	/** the reference to the problem */
	const Problem& __p;

	/** the coefficient matrix */
	Matrix __A;

	/** constant vector */
	Vector __b;

	/** iteration number */
	scalar __iter_n;

	/** rows of __A */
	Vector __Ar;

	/** elements of __b */
	scalar __be;

	/** parameter of Least Mean Square method */
	scalar __mu;

	/** the final result */
	Vector __x;
	//%} end of variables

	// void init();
	LeastMeanSquareSolver();

	/* "  something happens in the beginning of solving  " */
	void solvingBegin();

	/* "  the real solving part  " */
	void doSolve();

	/* "  the real computation of pre-treatment part  " */
	void doCompute();

	/* "    do one iteration
	 Returns:
	 The estimated error, a 'scalar'
	 " */
	scalar doOneIteration();

public:

//	typedef leastmeansquare_solver              ObjectCategory;
	/** constructor and deconstructor */
	//%{
	/* "    Initialization of the 'LeastMeanSquareSolver'.
	 Parameters:
	 s:    The least squares problem
	 mu:   The step size of least mean squares algorithm. 0.01 is used as default.
	 " */
	LeastMeanSquareSolver(const Problem& s, const scalar mu = 0.01);

	//%} end of constructor and deconstructor

	/** setter and getter */
	//%{
	/* "    Get the final result of 'LeastMeanSquareSolver'. 
	 Returns:
	 The final result of the solver, a 'Vector'
	 " */
	const Vector& result() const;
	//%}

	/* "    Calculate the square of l2-norm of 'Ax - b' in order to see how well the problem is solved using least mean squares algorithm.
	 Returns:
	 The square of l2-norm of 'Ax - b', a 'scalar'
	 " */
	scalar objective() const;

};

/*	"""     'LeastSquareSolver' is one of the solvers of least squares problems. 
 'LeastSquareSolver' solves the least squares problems by solving the equation 'A^TAx = A^Tb' directly.
 This class derives directly from another class called 'GeneralSolver' and takes 'Problem' and 'Time' as template.
 """ */

template<class Problem, class Time = NoTimeStatistics>
class LeastSquareSolver: public GeneralSolver<typename Problem::KernelTrait,
		Time> {
private:
	/*** the type of index used in 'LeastSquareSolver' ***/
	typedef typename Problem::KernelTrait::index index;
	/*** the type of scalar used in 'LeastSquareSolver' ***/
	typedef typename Problem::KernelTrait::scalar scalar;
	/*** the type of vector used in 'LeastSquareSolver' ***/
	typedef typename Problem::KernelTrait::Vector Vector;
	/*** the type of matrix used in 'LeastSquareSolver' ***/
	typedef typename Problem::KernelTrait::Matrix Matrix;

	/** private variables */
	//%{
	/** the reference to the problem */
	const Problem& __p;

	/** the coefficient matrix */
	Matrix __A;

	/** constant vector */
	Vector __b;

	/** the final result */
	Vector __x;

	//%} end of variables

	// void init();
	LeastSquareSolver();

	/* "  something happens in the beginning of solving  " */
	void solvingBegin();

	/* "  the real solving part  " */
	void doSolve();

	/* "  the real computation of pre-treatment part  " */
	void doCompute();

	/* "    do one iteration
	 Returns:
	 The estimated error, a 'scalar'
	 " */
	scalar doOneIteration();

public:

//	typedef leastsquare_solver              ObjectCategory;
	/** constructor and deconstructor */
	//%{
	/* "    Initialization of the 'LeastSquareSolver'.
	 Parameters:
	 s:    The least squares problem
	 " */
	LeastSquareSolver(const Problem& s);

	//%} end of constructor and deconstructor

	/** setter and getter */
	//%{
	/* "    Get the final result of 'LeastSquareSolver'. 
	 Returns:
	 The final result of the solver, a 'Vector'
	 " */
	const Vector& result() const;
	//%}

	/* "    Calculate the square of l2-norm of 'Ax - b' in order to see how well the problem is solved using 'LeastSquareSolver'.
	 Returns:
	 The square of l2-norm of 'Ax - b', a 'scalar'
	 " */
	scalar objective() const;

};

/*	"""     'RecursiveLeastSquareSolver' is one of the solvers of least squares problems using recursive least squares algorithm. 
 This class derives directly from another class called 'GeneralSolver' and takes 'Problem' and 'Time' as template.
 """ */

template<class Problem, class Time = NoTimeStatistics>
class RecursiveLeastSquareSolver: public GeneralSolver<
		typename Problem::KernelTrait, Time> {
private:
	/*** the type of index used in 'RecursiveLeastSquareSolver' ***/
	typedef typename Problem::KernelTrait::index index;
	/*** the type of scalar used in 'RecursiveLeastSquareSolver' ***/
	typedef typename Problem::KernelTrait::scalar scalar;
	/*** the type of vector used in 'RecursiveLeastSquareSolver' ***/
	typedef typename Problem::KernelTrait::Vector Vector;
	/*** the type of matrix used in 'RecursiveLeastSquareSolver' ***/
	typedef typename Problem::KernelTrait::Matrix Matrix;

	/** private variables */
	//%{
	/** the reference to the problem */
	const Problem& __p;

	/** the coefficient matrix */
	Matrix __A;

	/** constant vector */
	Vector __b;

	/** iteration number , also rows of A */
	scalar __iter_n;

	/** cols number of A */
	scalar __iter_m;

	/** rows of __A */
	Vector __Ar;

	/** elements of __b */
	scalar __be;

	/** parameter of Recursive Least Square method */
	scalar __lam;

	/** parameter of Recursive Least Square method */
	scalar __delta;

	/** intermediate variable */
	Matrix __inter_p;

	/** the final result */
	Vector __x;

	//%} end of variables

	// void init();
	RecursiveLeastSquareSolver();

	/* "  something happens in the beginning of solving  " */
	void solvingBegin();

	/* "  the real solving part  " */
	void doSolve();

	/* "  the real computation of pre-treatment part  " */
	void doCompute();

	/* "    do one iteration
	 Returns:
	 The estimated error, a 'scalar'
	 " */
	scalar doOneIteration();

public:

//	typedef recursiveleastsquare_solver              ObjectCategory;
	/** constructor and deconstructor */
	//%{
	/* "    Initialization of the 'RecursiveLeastSquareSolver'.
	 Parameters:
	 s:     The least squares problem
	 lam:   A parameter related to forgetting factor beta. 1 is used as default.
	 delta: A parameter used to initialize matrix P in recursive least squares algorithm. 250 is used as default.
	 " */
	RecursiveLeastSquareSolver(const Problem& s, const scalar lam = 1,
			const scalar delta = 250);

	//%} end of constructor and deconstructor

	/** setter and getter */
	//%{
	/* "    Get the final result of 'RecursiveLeastSquareSolver'. 
	 Returns:
	 The final result of the solver, a 'Vector'
	 " */
	const Vector& result() const;
	//%}

	/* "    Calculate the square of l2-norm of 'Ax - b' in order to see how well the problem is solved using recursive least squares algorithm.
	 Returns:
	 The square of l2-norm of 'Ax - b', a 'scalar'
	 " */
	scalar objective() const;

};

/*	"""     'LeastSquaresProblem' is one of the problem classes in open source library COPT. It describes a least square problem. 
 Briefly, a matrix 'A' and a target vector 'b' are given. We are supposed to find the optimal solution of 'Ax = b'.
 In other words, we need to minimize the l2-norm of 'Ax - b'.
 This class derives directly from another class called 'VectorProblem' and also takes 'kernal' as template.
 """ */
template<class kernel>
class LeastSquaresProblem: public VectorProblem<kernel> {
private:
	/*** the type of scalar used in 'LeastSquaresProblem' ***/
	typedef typename kernel::scalar scalar;
	/*** the type of vector used in 'LeastSquaresProblem' ***/
	typedef typename kernel::Vector Vector;
	/*** the type of matrix used in 'LeastSquaresProblem' ***/
	typedef typename kernel::Matrix Matrix;

	/** the coefficient matrix */
	const Matrix& __A;
	/** the target vector */
	const Vector& __b;

public:
	/*** the kernel that is used in the solver ***/
	typedef kernel KernelTrait;

	// typedef leastsquares_problem					ObjectCategory;

	/* "    Initialization of the least squares problem.
	 Parameters:
	 A:    The given coefficient matrix
	 b:    The given target vector
	 " */
	LeastSquaresProblem(const Matrix& A, const Vector& b);

	// void leastMeanSquareSolve();
	// void leastSquareSolve();
	// void recursiveLeastSquareSolve();

	/* "    Check whether the input is valid. This means that the column number of the given coefficient matrix 'A' must be equal to the size of the parameter vector 'x' to be find.
	 Parameters:
	 x:    The parameter vector 'x' to be find
	 Returns:
	 True if the input is valid, false otherwise.
	 " */
	bool isValidInput(const Vector& x) const;

	/* "    Check whether the problem is a valid problem. This means that the row number of the given coefficient matrix 'A' must be equal to the size of the given target vector 'b'.
	 Returns:
	 True if the problem is a valid one, false otherwise.
	 " */
	bool isValid() const;

	/** getter and setter */
	//%{
	/* "    Get the coefficient matrix 'A' of a least squares problem.
	 Returns:
	 The coefficient matrix A
	 " */
	const Matrix& matA() const;

	/* "    Get the target vector 'b' of a least squares problem.
	 Returns:
	 The target vector 'b'
	 " */
	const Vector& obB() const;

	/* "    Calculate the square of l2-norm of 'Ax - b' in order to see how well the problem is solved.
	 Parameters:
	 x:    The parameter vector 'x'
	 Returns:
	 The square of l2-norm of 'Ax - b', a 'scalar'
	 " */
	scalar objective(const Vector& x) const;
	//%}
};

/*********************Implementation of LeastMeanSquareSolver ******************/

template<class Problem, class Time>
void LeastMeanSquareSolver<Problem, Time>::doSolve() {
	this->__iter_num = 0;
	do {
		this->__estimated_error = this->oneIteration();
		if (++this->__iter_num >= __iter_n) {
			this->__terminal_type = this->MaxIteration;
			break;
		}
	} while (this->__estimated_error > this->__thresh);
	if (this->__estimated_error <= this->__thresh)
		this->__terminal_type = this->Optimal;
}

template<class Problem, class Time>
void LeastMeanSquareSolver<Problem, Time>::doCompute() {
	__A = __p.matA();
	__b = __p.obB();
}

template<class Problem, class Time>
LeastMeanSquareSolver<Problem, Time>::LeastMeanSquareSolver(const Problem& s,
		const scalar mu) :
		__p(s), __mu(mu) {
	this->doCompute();
}

template<class Problem, class Time>
void LeastMeanSquareSolver<Problem, Time>::solvingBegin() {
	__x.resize(__A.cols());
	__iter_n = __A.rows();
}

template<class Problem, class Time>
typename LeastMeanSquareSolver<Problem, Time>::scalar LeastMeanSquareSolver<
		Problem, Time>::doOneIteration() {
	Vector xp = __x;
	__Ar = __A.row(this->__iter_num);
	__be = __b[this->__iter_num];
	scalar e = __be - __Ar.dot(__x);
	__x = __x + __mu * e * __Ar;
	return (__A * __x - __b).squaredNorm();
}

template<class Problem, class Time>
const typename LeastMeanSquareSolver<Problem, Time>::Vector& LeastMeanSquareSolver<
		Problem, Time>::result() const {
	return __x;
}

template<class Problem, class Time>
typename LeastMeanSquareSolver<Problem, Time>::scalar LeastMeanSquareSolver<
		Problem, Time>::objective() const {
	return __p.objective(__x);
}

/*********************Implementation of LeastSquareSolver ******************/

template<class Problem, class Time>
void LeastSquareSolver<Problem, Time>::doSolve() {
	this->__estimated_error = this->oneIteration();
	if (this->__estimated_error <= this->__thresh)
		this->__terminal_type = this->Optimal;
}

template<class Problem, class Time>
void LeastSquareSolver<Problem, Time>::doCompute() {
	__A = __p.matA();
	__b = __p.obB();
}

template<class Problem, class Time>
LeastSquareSolver<Problem, Time>::LeastSquareSolver(const Problem& s) :
		__p(s) {
	this->doCompute();
}

template<class Problem, class Time>
void LeastSquareSolver<Problem, Time>::solvingBegin() {
	__x.resize(__A.cols());
}

template<class Problem, class Time>
typename LeastSquareSolver<Problem, Time>::scalar LeastSquareSolver<Problem,
		Time>::doOneIteration() {
	Vector xp = __x;
	__x = (__A.transpose() * __A).solve(__A.transpose() * __b);
	return (__A * __x - __b).squaredNorm();
}

template<class Problem, class Time>
const typename LeastSquareSolver<Problem, Time>::Vector& LeastSquareSolver<
		Problem, Time>::result() const {
	return __x;
}

template<class Problem, class Time>
typename LeastSquareSolver<Problem, Time>::scalar LeastSquareSolver<Problem,
		Time>::objective() const {
	return __p.objective(__x);
}

/*********************Implementation of RecursiveLeastSquareSolver ******************/

template<class Problem, class Time>
void RecursiveLeastSquareSolver<Problem, Time>::doSolve() {
	this->__iter_num = 0;
	do {
		this->__estimated_error = this->oneIteration();
		if (++this->__iter_num >= __iter_n) {
			this->__terminal_type = this->MaxIteration;
			break;
		}
	} while (this->__estimated_error > this->__thresh);
	if (this->__estimated_error <= this->__thresh)
		this->__terminal_type = this->Optimal;
}

template<class Problem, class Time>
void RecursiveLeastSquareSolver<Problem, Time>::doCompute() {
	__A = __p.matA();
	__b = __p.obB();
}

template<class Problem, class Time>
RecursiveLeastSquareSolver<Problem, Time>::RecursiveLeastSquareSolver(
		const Problem& s, const scalar lam, const scalar delta) :
		__p(s), __lam(lam), __delta(delta) {
	this->doCompute();
}

template<class Problem, class Time>
void RecursiveLeastSquareSolver<Problem, Time>::solvingBegin() {
	__x.resize(__A.cols());
	__iter_n = __A.rows();
	__iter_m = __A.cols();
	__inter_p = __delta * Matrix::identity(__iter_m, __iter_m);
}

template<class Problem, class Time>
typename RecursiveLeastSquareSolver<Problem, Time>::scalar RecursiveLeastSquareSolver<
		Problem, Time>::doOneIteration() {
	Vector xp = __x;
	__Ar = __A.row(this->__iter_num);
	__be = __b[this->__iter_num];
	Vector pai = __inter_p.transpose() * __Ar;
	scalar gama = __lam + pai.dot(__Ar);
	Vector k = pai * (1.0 / gama);
	scalar alpha = __be - __x.dot(__Ar);
	__x = __x + k * alpha;
	Matrix pp = k.mulTrans(pai);
	__inter_p = (1.0 / __lam) * (__inter_p - pp);
	return (__A * __x - __b).squaredNorm();
}

template<class Problem, class Time>
const typename RecursiveLeastSquareSolver<Problem, Time>::Vector& RecursiveLeastSquareSolver<
		Problem, Time>::result() const {
	return __x;
}

template<class Problem, class Time>
typename RecursiveLeastSquareSolver<Problem, Time>::scalar RecursiveLeastSquareSolver<
		Problem, Time>::objective() const {
	return __p.objective(__x);
}

/***********************Implementation of LeastSquaresProblem************************/

template<class kernel>
LeastSquaresProblem<kernel>::LeastSquaresProblem(const Matrix& A,
		const Vector& b) :
		__A(A), __b(b) {
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
bool LeastSquaresProblem<kernel>::isValid() const {
	if (__A.rows() != __b.size())
		return false;
	else
		return true;
}

template<class kernel>
bool LeastSquaresProblem<kernel>::isValidInput(const Vector& x) const {
	if (__A.cols() != x.size())
		return false;
	else
		return true;
}

template<class kernel>
const typename kernel::Matrix& LeastSquaresProblem<kernel>::matA() const {
	return __A;
}

template<class kernel>
const typename kernel::Vector& LeastSquaresProblem<kernel>::obB() const {
	return __b;
}

template<class kernel>
typename kernel::scalar LeastSquaresProblem<kernel>::objective(
		const Vector& x) const {
	return (__A * x - __b).squaredNorm();
}
}

#endif
