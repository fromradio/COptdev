//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU
#ifndef BFGS_H
#define BFGS_H

namespace COPT{
template<class VFunc>
void BFGSMethod(
	const VFunc& 						func,
	const typename VFunc::ScalarType 	c1,
	const typename VFunc::ScalarType 	c2,
	const typename VFunc::ScalarType 	sigma,
	typename VFunc::Vector& 			x,
	typename VFunc::ScalarType& 		tol_error,
	int& 								iters
	)
{
	typedef typename VFunc::ScalarType 		Scalar;
	typedef typename VFunc::Vector 			Vector;
	typedef typename VFunc::Matrix 			Matrix;

	// gradient;
	Vector gradient = func.gradient(x);
	// Vector xformer = x;
	Vector gradientformer = gradient;

	Scalar tol = tol_error*tol_error;
	tol_error = gradient.squaredNorm();
	int maxIter = iters;

	Matrix H = Matrix::identity(x.size(),x.size(),std::sqrt(tol_error)*sigma);
	Matrix I = Matrix::identity(x.size(),x.size());
	while (tol_error>tol){
		int numbers = 100;
		Scalar steplength = 1.0;
		Vector direction = -(H*gradient);
		backTrackingWithWolfeCondition(func,x,gradient,direction,1.0,c1,c2,steplength,numbers);
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
		if(iters>=maxIter)
			break;
	}
}
}// End of namespace COPT

#endif