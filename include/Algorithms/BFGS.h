//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU
#ifndef BFGS_H
#define BFGS_H

template<class VFunc>
void BFGSMethod(
	const VFunc& 						func,
	const typename VFunc::ScalarType 	c1,
	const typename VFunc::ScalarType 	c2,
	const typename VFunc::ScalarType 	sigma,
	typename VFunc::Vector& 			x,
	typename VFunc::ScalarType& 		tol_error,
	int& 								iters,
	)
{
	typedef typename VFunc::ScalarType 		Scalar;
	typedef typename VFunc::Vector 			Vector;
	typedef typename VFunc::Matrix 			Matrix;

	// gradient;
	Vector gradient = func.gradient(x);
	
}

#endif