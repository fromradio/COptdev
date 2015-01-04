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


#ifndef CONJUGATE_GRADIENT_METHOD_HPP__
#define CONJUGATE_GRADIENT_METHOD_HPP__




namespace COPT{
/*
 *		A function solving linear conjugate gradient problem
 *		/param mat, the left hand matrix
 *		/param rhs, the right hand side vector
 *		/param x, the initial on input and solution on output
 *		/param iters, the maximum iteration on input and the real iteration on output
 *		/param tol, the tolerance error on input and estimated error on output
 *
 *		Notice:
 *			estimated error is computed as residual.squareNorm()/rhs.squaredNorm()
 */
template <class Matrix,class Vector>
void conjugateGradientWithoutPrecondition(
	const Matrix& mat,
	const Vector& rhs,
	Vector& x, 
	int& iters,
	typename Vector::ScalarType& tol_error)
{
	using std::sqrt;
	typedef typename Vector::ScalarType			ScalarType;

	int maxIter = iters;

	Vector residual(mat*x-rhs);
	Vector p(-residual);

	ScalarType residualnorm2 = residual.squaredNorm();
	ScalarType rhsnorm2 	 = rhs.squaredNorm();
	ScalarType tol = tol_error*rhsnorm2;
	iters = 0;
	while (residualnorm2 > tol) {
		ScalarType alpha = residualnorm2/(p.dot(mat*p));
		x = x + alpha*p;
		residual = residual + alpha*(mat*p);
		ScalarType formalnorm = residualnorm2;
		residualnorm2 = residual.squaredNorm();
		ScalarType beta = residualnorm2/formalnorm;
		p = -residual + beta*p;
		++ iters;
		if (iters>=maxIter)
			break;
	}

	tol_error = sqrt(residualnorm2/rhsnorm2);
}

/*			conjugate gradient method for non-linear problem
 */
// template<class Function class Vector>
// void nonlinearConjugateGradient(
// 	const Function& func,
// 	Vector& x)
// {

// }
} // End of namespace COPT



#endif