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

#ifndef NON_LINEAR_SOLVER_HPP__
#define NON_LINEAR_SOLVER_HPP__

namespace COPT {
/*			Numerical solver for non-linear problems
 *			Currently, three methods are contained:
 *			Steepest descent method
 *			Newton's method
 *			BFGS method
 */
template<class VFunc>
class NonLinearSolver {
public:
	enum SolverType {
		SDM,		// steepest descent method
		NM,			// Newton's method
		BFGS		// BFGS method
	};
private:
	typedef typename VFunc::ScalarType ScalarType;
	typedef typename VFunc::Vector Vector;
	// 	the reference to target function
	const VFunc& __func;
	//	the tolerance of error
	ScalarType __tol;
	//	the maximum number of iterations
	int __maxiternum;
	//	the first constant for back tracking method
	ScalarType __c1;
	//	the second constant for wolfe condition judging (back tracking method)
	ScalarType __c2;
	//	the ratio for back tracking method
	ScalarType __ratio;
	//	the sigma value for BFGS method
	ScalarType __sigma;
	//	the maximum number for back tracking
	int __tracknum;
	//	the final error
	ScalarType __error;
	//	the final number of iterations
	int __iters;
	//	the final result
	Vector __x;
	//	the status of solver, true for converged and false for un-converged
	bool __status;
	//	the type of the solver, SDM as default
	SolverType __type;

	/*
	 private functions
	 */
	//	the main algorithm
	bool doSolve();
public:
	NonLinearSolver(const VFunc& func, const ScalarType tol = 1e-5,
			const ScalarType maxiternum = 1000, const ScalarType c1 = 1e-3,
			const ScalarType c2 = 0.5, const ScalarType ratio = 0.7,
			const ScalarType sigma = 1.0, const int __tracknum = 100,
			const SolverType type = SDM);

	//		set the type of the solver: SDM, NM or BFGS
	void setType(const SolverType type);

	//		set the maximus iteration number
	void setIterationNum(const int maxiternum);

	//		set the error threshold
	void setErrorThreshold(const ScalarType tol);

	//		interface for solving the problem
	void solve(const Vector& vec);

	//		print the solver information
	void printInfo();
};
}

#endif
