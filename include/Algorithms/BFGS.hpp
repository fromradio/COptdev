//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU
#ifndef BFGS_H
#define BFGS_H

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