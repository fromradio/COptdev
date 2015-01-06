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


#ifndef LINE_SEARCH_HPP__
#define LINE_SEARCH_HPP__


namespace COPT
{

/*
 *		judge whether wolfe condition is satisfied
 *		first condition is that
 *			f(x+alpha*p)<=f(x)+c1*alpha*\nablaf(x).dot(p)
 *		second condition is that
 *			\nablaf(x+alpha*p).dot(p)>=c2*\nablaf(x).dot(p)
 *		/param func:			the function
 *		/param x:				current point
 *		/param alpha:			the step length
 *		/param p:				the direction
 *		/param gradient:		the direction of gradient
 *		/param c1:				the first constant
 *		/param c2:				the second constant
 */
template<class VFunc>
bool judgeWolfeCondition(
	const VFunc& func,
	const typename VFunc::Vector& x,
	const typename VFunc::ScalarType alpha,
	const typename VFunc::Vector& p,
	const typename VFunc::Vector& gradient,
	const typename VFunc::ScalarType c1,
	const typename VFunc::ScalarType c2 )
{
	typename VFunc::ScalarType dot = gradient.dot(p);
	if(func(x+alpha*p)<=func(x)+c1*alpha*dot)
		if(func.gradient(x+alpha*p).dot(p)>=c2*dot)
			return true;
		else
			return false;
	else
		return false;
}

/*
 *		judge whether strong Wolfe condition is satisfied
 *		first condition is that
 *			f(x+alpha*p)<=f(x)+c1*alpha*\nablaf(x).dot(p)
 *		second condition is that
 *			abs(\nablaf(x+alpha*p).dot(p))>=abs(c2*\nablaf(x).dot(p))
 *		/param func:			the function
 *		/param x:				current point
 *		/param alpha:			the step length
 *		/param p:				the direction
 *		/param gradient:		the direction of gradient
 *		/param c1:				the first constant
 *		/param c2:				the second constant
 */
template<class VFunc>
bool judgeStrongWolfeCondition(
	const VFunc& func,
	const typename VFunc::Vector& x,
	const typename VFunc::ScalarType alpha,
	const typename VFunc::Vector& p,
	const typename VFunc::Vector& gradient,
	const typename VFunc::ScalarType c1,
	const typename VFunc::ScalarType c2)
{
	typename VFunc::ScalarType dot = gradient.dot(p);
	if(func(x+alpha*p)<=func(x)+c1*alpha*dot)
		if(std::abs(func.gradient(x+alpha*p).dot(p))<=c2*std::abs(dot))
			return true;
		else
			return false;
	else
		return false;
}


/*
 *		Backtracking method to inexactly find the step length
 *		/param func 		the given function
 *		/param x 			current iteration point
 *		/param gradient		the gradient of function at current point
 *		/param direction	the descent direction at the point
 *		/param rho			the shrink ratio
 *		/param c 			the constant ratio
 *		/param alpha		the initial alpha on input and estimated step length on output
 *		/param iters 		the max iteration number on input and final iteration number on output
 */
template<class Function,class Vector>
void findStepLengthBackTracking(
	const Function& func,
	const Vector& x,
	const Vector& gradient,
	const Vector& direction,
	typename Vector::ScalarType rho,
	typename Vector::ScalarType c,
	typename Vector::ScalarType& alpha,
	int& iters)
{
	typedef typename Vector::ScalarType 		ScalarType;

	int maxIter = iters;
	iters = 0;
	ScalarType f = func(x);
	ScalarType dot = gradient.dot(direction);
	std::cout<<gradient<<std::endl;
	std::cout<<dot<<std::endl;
	while(iters < maxIter){
		if(func(x+alpha*direction)<=f+c*alpha*dot)
			break;
		else
			alpha = rho*alpha;
		++ iters;
	}
	std::cout<<iters<<' '<<alpha<<std::endl;
	std::cout<<x+alpha*direction<<std::endl;
}


/*		Back tracking method to find the step length satisfying Wolfe condition
 *		/param func:				the input function
 *		/param x:					current point
 *		/param gradient:			the gradient of current point
 *		/param direction:			the descent direction
 *		/param rho:					the scaling ratio
 *		/param c1:					the first constant
 *		/param c2:					the second constant
 *		/param alpha:				the initial step length on input and result on output
 */
template<class Function,class Vector>
void backTrackingWithWolfeCondition(
	const Function& func,
	const Vector& x,
	const Vector& gradient,
	const Vector& direction,
	const typename Vector::ScalarType rho,
	const typename Vector::ScalarType c1,
	const typename Vector::ScalarType c2,
	typename Vector::ScalarType& alpha,
	int& iters)
{
	typedef typename Vector::ScalarType 		ScalarType;

	int maxIter = iters;
	iters = 0;
	ScalarType f = func(x);
	ScalarType dot = gradient.dot(direction);
	while(iters < maxIter){
		if(func(x+alpha*direction)<=f+c1*alpha*dot)
			if(func.gradient(x+alpha*direction).dot(direction)>=c2*dot)
				break;
			else
				;
		else
			alpha = rho*alpha;
		++ iters;
	}
}


/*
 *		Steepest Descent approach solving non-linear problem
 *		/param func 		input function
 *		/param rho 			ratio for finding step length
 *		/param c 			the constant ration
 *		/param x 			the inital x on input and optimized x on output
 *		/param tol_error	the tolerance on input and estimated error on output
 *		/param iters		maximum iteration number on input and real iterations on output 
 */
template<class Function>
void steepestDescentUsingBackTracking( 
	const Function& func , 
	const typename Function::ScalarType rho, 
	const typename Function::ScalarType c , 
	typename Function::Vector& x ,
	typename Function::ScalarType& tol_error, 
	int& iters ,
	const int tracknum = 100
	)
{
	using std::sqrt;
	typedef typename Function::ScalarType 			ScalarType;
	typedef typename Function::Vector 				Vector;

	int maxIter = iters;
	iters = 0;
	ScalarType tol = tol_error*tol_error;
	Vector gradient = func.gradient(x);
	Vector direction = -gradient;
	ScalarType error = direction.squaredNorm();
	while (error>tol){
		ScalarType steplength = 1;
		int biter = tracknum;
		findStepLengthBackTracking(func,x,gradient,direction,rho,c,steplength,biter);
		// backTrackingWithWolfeCondition(func,x,gradient,direction,rho,c,0.4,steplength,biter);
		x = x + steplength*direction;
		++ iters;
		if(iters>=maxIter)
			break;
		gradient = func.gradient(x);
		direction = -gradient;
		error = gradient.squaredNorm();
	}

	tol_error = sqrt(error);
}


/*		Newton method
 *		/param func:			the input function
 *		/param x:				intial point on input and result on output
 *		/param tol_error:		tolerance on input and final estimated error on output
 *		/param iters:			maximum iteration number on input and final iteration number on output
 */
template<class VFunc>
void newtonMethod(
	const VFunc& func,
	typename VFunc::Vector& x,
	typename VFunc::ScalarType& tol_error,
	int & iters)
{
	typedef typename VFunc::Vector		Vector;
	typedef typename VFunc::ScalarType	ScalarType;
	typedef typename VFunc::Matrix 		Matrix;
	int maxIter = iters;
	iters = 0;
	Vector gradient = func.gradient(x);
	Vector direction;
	ScalarType tol = tol_error*tol_error;
	tol_error = gradient.squaredNorm();
	while(tol_error>tol){
		Matrix hessian = func.hessian(x);
		direction = hessian.solve(gradient);
		x = x - direction;
		gradient = func.gradient(x);
		++ iters;
		tol_error = gradient.squaredNorm();
		if ( iters >=  maxIter)
			break;
	}
	tol_error = std::sqrt(tol_error);
}

}// End of namespace COPT

#endif