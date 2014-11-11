// 		This file is part of open library COPT
//		Copyright (c) MathU
//		Written by Ruimin Wang, ruimin.wang13@gmail.com

#ifndef LINE_SEARCH_H
#define LINE_SEARCH_H


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
		if(func.gradient(x+alpha*p).dot(p)>=dot)
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
	while(iters < maxIter){
		if(func(x+alpha*direction)<=f+c*alpha*dot)
			break;
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
template<class Function,class Vector>
void steepestDescentUsingBackTracking( 
	const Function& func , 
	const typename Vector::ScalarType rho, 
	const typename Vector::ScalarType c , 
	Vector& x ,
	typename Vector::ScalarType& tol_error, 
	int& iters )
{
	using std::sqrt;
	typedef typename Vector::ScalarType 			ScalarType;

	int maxIter = iters;
	iters = 0;
	ScalarType tol = tol_error*tol_error;
	Vector gradient = func.gradient(x);
	Vector direction = -gradient;
	ScalarType error = direction.squaredNorm();
	while (error>tol){
		ScalarType steplength = 1;
		int biter = 100;
		findStepLengthBackTracking(func,x,gradient,direction,rho,c,steplength,biter);
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
	ScalarType error = gradient.squaredNorm();
	while(error>tol){
		Matrix hessian = func.hessian(x);
		direction = hessian.solve(gradient);
		x = x - direction;
		gradient = func.gradient(x);
		++ iters;
		error = gradient.squaredNorm();
		if ( iters >=  maxIter)
			break;
	}
}

template<class VFunc>
void BFGSMethod(
	)
{
	
}
};

#endif