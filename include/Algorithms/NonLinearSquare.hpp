// This file is part of COPT, a light-weight C++ based optimization open source library
//
// Copyright (C) 2015 Ruimin Wang <ruimin.wang13@gmail.com>
// Copyright (C) 2015 MathU
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

#ifndef NON_LINEAR_SQUARE_HPP__
#define NON_LINEAR_SQUARE_HPP__

namespace COPT {
/* 	This algorithm is used to solve non-linear least squares problems.
 *	It is a uncontrained optimization problem. Here, we use Levenberg
 *	Marquardt Algorithm(LMA) to solve the problem. The LMA interpolates 
 *	between the Gauss–Newton algorithm (GNA) and the method of
 *	gradient descent. The LMA is more robust than the GNA. 
 *	LMA can also be viewed as Gauss–Newton using a trust region approach.
 */

/*		class of non-linear square algorithm
 *
 */
template<class Problem, class Time = NoTimeStatistics>
class NonLinearSquare: public GeneralSolver<typename Problem::KernelTrait, Time> {
private:
	typedef typename Problem::KernelTrait::index index;
	typedef typename Problem::KernelTrait::scalar scalar;
	typedef typename Problem::KernelTrait::Vector Vector;
	typedef typename Problem::KernelTrait::Matrix Matrix;

	/*
	 *		vector functions that are used
	 *				
	 */

	const Problem& __p;
	//	point value
	Vector __x;
	//	function dimension
	index __m;
	//	variable dimension
	index __n;
	//	constants
	scalar __tao;
	scalar __v;
	//	iteration
	index __maxiteration;
	/*	A is approximately equal to Hessian	matrix 	*/
	Matrix __A;
	/*	g is approximately equal to Gradient	*/
	Vector __g;
	//	scaling factor
	scalar __u;
	bool __found;

	/*
	 *		non-linear square solver
	 */

	void solvingBegin();
	void doCompute();
	scalar doOneIteration();
	void doSolve();
public:
	/*
	 *		Constructor
	 */
	NonLinearSquare(const Problem& p, const Vector& x, const scalar tao = 1e-3,
			const scalar v = 2.0, const index iternum = 1000);

	// 	output
	void printInfo();
	//	objective value
	scalar objective() const;
	// 	optimal point
	const Vector& result() const;

};

/***********************Implementation of NonLinearSquare (LM)************************/

template<class Problem, class Time>
void NonLinearSquare<Problem, Time>::doCompute() {
	__A = __p.jacobi(__x).transpose() * (__p.jacobi(__x));
	__g = __p.jacobi(__x).transpose() * (__p.value(__x));

	scalar max = fabs(__A(0, 0));
	for (int i = 0; i < __n; i++) {
		if (fabs(__A(i, i)) > max)
			max = __A(i, i);
	}
	__u = __tao * max;
}

template<class Problem, class Time>
NonLinearSquare<Problem, Time>::NonLinearSquare(const Problem& p,
		const Vector& x, const scalar tao, const scalar v, const index iternum) :
		__p(p), __x(x), __tao(tao), __v(v), __maxiteration(iternum) {

	__m = __p.functionDimension();
	__n = __p.variableDimension();
	this->doCompute();
}

template<class Problem, class Time>
void NonLinearSquare<Problem, Time>::solvingBegin() {
	__m = __p.functionDimension();
	__n = __p.variableDimension();
}

template<class Problem, class Time>
typename NonLinearSquare<Problem, Time>::scalar NonLinearSquare<Problem, Time>::doOneIteration() {
	/*
	 Levenberg Marquardt method
	 */

	scalar e = 1e-8;
	scalar rho;
	Vector h;
	Vector x;

	__found = (sqrt(__g.dot(__g)) <= e);
	Matrix B(__A);
	Vector k(__g);

	for (int i = 0; i < __n; i++) {
		B(i, i) = __A(i, i) + __u;
		k[i] = -__g[i];
	}

	/*	Solve linear equations	*/
	h = B.solve(k);

	if (h.dot(h) < e * (e + __x.dot(__x))) {
		__found = true;
		return sqrt(h.dot(h));
	} else {
		x = __x + h;

		/*	Judgment flag	rho = (F(__x)-F(x))/L(h,u,g)	*/
		rho = (__p.value(__x).dot(__p.value(__x)) / 2
				- __p.value(x).dot(__p.value(x)) / 2)
				/ (0.5 * h.dot(__u * h - __g));

		if (rho > 0.0) {
			__x = x;
			__A = __p.jacobi(__x).transpose() * (__p.jacobi(__x));
			__g = __p.jacobi(__x).transpose() * (__p.value(__x));
			__found = (sqrt(__g.dot(__g)) <= e);
			scalar uu = 1 / 3;
			if (1 - (2 * rho - 1) * (2 * rho - 1) * (2 * rho - 1) > uu) {
				uu = 1 - (2 * rho - 1) * (2 * rho - 1) * (2 * rho - 1);
			}
			__u = __u * uu;
			__v = 2.0;
		} else {
			__u = __u * __v;
			__v = 2 * __v;
		}
		return sqrt(__g.dot(__g));
	}
}

template<class Problem, class Time>
void NonLinearSquare<Problem, Time>::doSolve() {
	this->__iter_num = 0;
	do {
		this->__estimated_error = this->oneIteration();
		if (++this->__iter_num >= __maxiteration) {
			this->__terminal_type = this->MaxIteration;
			break;
		}
	} while (!__found);
	if (__found)
		this->__terminal_type = this->Optimal;
}

template<class Problem, class Time>
void NonLinearSquare<Problem, Time>::printInfo() {
	std::cout << "The optimal point is " << __x << std::endl;
}

template<class Problem, class Time>
typename NonLinearSquare<Problem, Time>::scalar NonLinearSquare<Problem, Time>::objective() const {
	return __p.objective(__x);
}

template<class Problem, class Time>
const typename NonLinearSquare<Problem, Time>::Vector& NonLinearSquare<Problem,
		Time>::result() const {
	return __x;
}

/*		class of non-linear square problem
 *		min 0.5*|f(x)|^2 	There are no contrains
 *		f(x) = ( f1(x1, ... , xn) , ... , fm(x1, ... , xn) ) is a non-linear functin
 */
template<class kernel>
class NonLinearSquareProblem: public VectorProblem<kernel> {
private:
	typedef typename kernel::index index;
	typedef typename kernel::scalar scalar;
	typedef typename kernel::Matrix Matrix;
	typedef typename kernel::Vector Vector;

	VectorFunctionSystem<scalar> __vfs;
	index __m;
	index __n;

public:
	NonLinearSquareProblem(const VectorFunctionSystem<scalar>& vfs);

	bool isValidInput(const Vector& x) const;
	bool isValid() const;
	const index& functionDimension() const;
	const index& variableDimension() const;
	Vector value(const Vector& x) const;
	Matrix jacobi(const Vector& x) const;
	scalar objective(const Vector& x) const;

};

/***********************Implementation of NonLinearSquareProblem************************/

template<class kernel>
NonLinearSquareProblem<kernel>::NonLinearSquareProblem(
		const VectorFunctionSystem<scalar>& vfs) :
		__vfs(vfs) {
	__m = __vfs.functionDimension();
	__n = __vfs.variableDimension();
}

template<class kernel>
bool NonLinearSquareProblem<kernel>::isValidInput(const Vector& x) const {
	return (x.size() == __n);
}

template<class kernel>
bool NonLinearSquareProblem<kernel>::isValid() const {
	return (__m <= __n);
}

template<class kernel>
const typename kernel::index& NonLinearSquareProblem<kernel>::functionDimension() const {
	return __m;
}

template<class kernel>
const typename kernel::index& NonLinearSquareProblem<kernel>::variableDimension() const {
	return __n;
}

template<class kernel>
typename kernel::Vector NonLinearSquareProblem<kernel>::value(
		const Vector& x) const {
	return __vfs.functionValue(x);
}

template<class kernel>
typename kernel::Matrix NonLinearSquareProblem<kernel>::jacobi(
		const Vector& x) const {
	return __vfs.jacobiFunction(x);
}

template<class kernel>
typename kernel::scalar NonLinearSquareProblem<kernel>::objective(
		const Vector& x) const {
	return 0.5 * __vfs.functionValue(x).dot(__vfs.functionValue(x));
}

} // End of namespace COPT

#endif
