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


#ifndef BFGS_HPP__
#define BFGS_HPP__

namespace COPT{


/*				classical BFGS approach
 *				/param func:			input function
 *				/param c1:				c1 for wolfe condition
 *				/param c2:				c2 for wolfe condition
 *				/param sigma:			initial constant for the first H matrix
 *				/param x:				inital point on input and result on output
 *				/param tol_error:		tolerence on input and final error on output
 *				/param iters:			maximum iterations on input and final iteration number on output
 *				/param rho:				scaling ratio for Wolfe condition
 *				/param tracknum:		the maximum number using back tracking method to find the step length
 */
template<class VFunc>
void BFGSMethod(
	const VFunc& 						func,
	const typename VFunc::ScalarType 	c1,
	const typename VFunc::ScalarType 	c2,
	const typename VFunc::ScalarType 	sigma,
	typename VFunc::Vector& 			x,
	typename VFunc::ScalarType& 		tol_error,
	int& 								iters,
	const typename VFunc::ScalarType 	rho = 0.7,
	const int 							tracknum = 100
	)
{
	typedef typename VFunc::ScalarType 		Scalar;
	typedef typename VFunc::Vector 			Vector;
	typedef typename VFunc::Matrix 			Matrix;

	// gradient;
	Vector gradient = func.gradient(x);
	Vector gradientformer = gradient;

	Scalar tol = tol_error*tol_error;
	tol_error = gradient.squaredNorm();
	int maxIter = iters;
	iters = 0;

	Matrix H = Matrix::identity(x.size(),x.size(),std::sqrt(tol_error)*sigma);
	Matrix I = Matrix::identity(x.size(),x.size());
	while (tol_error>tol){
		int numbers = tracknum;
		Scalar steplength = 1.0;
		Vector direction = -(H*gradient);
		backTrackingWithWolfeCondition(func,x,gradient,direction,rho,c1,c2,steplength,numbers);
		Vector s = steplength*direction;
		x = x+s;
		gradientformer = gradient;
		gradient = func.gradient(x);
		Vector y = gradient-gradientformer;
		Scalar rho = 1.0/(y.dot(s));
		s = rho*s;
		H = (I-s.mulTrans(y))*H*(I-y.mulTrans(s))+1.0/rho*(s.mulTrans(s));
		tol_error = gradient.squaredNorm();
		++ iters;
		if(iters>=maxIter){
			tol_error = std::sqrt(tol_error);
			break;
		}
	}
	tol_error = std::sqrt(tol_error);
}
}// End of namespace COPT

#endif