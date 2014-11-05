// 		This file is part of open library COPT
//		Copyright (c) MathU
//		Written by Ruimin Wang, ruimin.wang13@gmail.com

#ifndef LINE_SEARCH_H
#define LINE_SEARCH_H


namespace COPT
{

/*
 *		Backtracking method to inexactly find the step length
 *		/params func 		the given function
 *		/params x 			current iteration point
 *		/params gradient	the gradient of function at current point
 *		/params direction	the descent direction at the point
 *		/params rho			the shrink ratio
 *		/params c 			the constant ratio
 *		/params alpha		the initial alpha on input and estimated step length on output
 *		/params iters 		the max iteration number on input and final iteration number on output
 */
template<class Function,class Vector>
void findStepLengthBackTracking(const Function& func,const Vector& x,const Vector& gradient,const Vector& direction,typename Vector::ScalarType rho,typename Vector::ScalarType c,typename Vector::ScalarType& alpha,int& iters)
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


template<class Function,class Vector>
void steepestDescentUsingBackTracking( const Function& func , typename Vector::ScalarType rho, typename Vector::ScalarType c , Vector& x ,typename Vector::ScalarType& tol_error, int& iters )
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
};

#endif