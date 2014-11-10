//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

// Function.h
#ifndef FUNCTION_H
#define FUNCTION_H
/*
	Second version of 'Functions' classes
*/
namespace COPT{

/*
 *	Scalar function
 *		input: scalar x
 *		output: scalar y
 *	template
 *		T stands for the float type that is used, for example double or float
 *
 */
template<class T>
class ScalarFunction{
protected:

public:
	typedef 		T 				ScalarType;
	typedef 		T 			 	FT;
	ScalarFunction (  ){}
	// the deconstructor
	virtual ~ScalarFunction() {}
	// basic operation computing the value of the function
	virtual FT operator() 	( FT x ) const = 0;
	// basic operation computing the differential of the function
	virtual FT diff  		( FT x ) const {
		return ScalarDifferential<ScalarFunction>(*this).diff(x);
	}
	virtual FT sDiff 		( FT x ) const {
		return ScalarDifferential<ScalarFunction>(*this).sDiff(x);
	}
};

/*
 *	'VectorFunction' returns a scalar value of a function defined on a vector
 *	call statement:
 *		double value = VectorFunction(vec):
 *			\param 'vec': input type Vector
 *			Returns the result of the function
 *		Vector diff = VectorFunction.gradient(vec):
 *			\param 'vec': input type Vector
 *			Returns the gradient of the function
 *		Matrix hessian = VectorFunction.hessian(vec):
 *			\param 'vec': input type Vector
 *			Returns the Hessian matrix of the function
 */
template<class VT>
class VectorFunction{
protected:
	// the dimension of the problem
	int 		__dim;
public:
	typedef 			VT 								Vector;
	typedef typename 	Vector::ScalarType 				ScalarType;
	typedef 			Matrix<ScalarType>				Matrix;
	typedef typename 	Vector::ScalarType 				FT;

	VectorFunction ( ):__dim(0) {}

	virtual ~VectorFunction() {}

	int 				dimension () {return __dim;}

	virtual FT 			operator() (const Vector& vec ) const = 0;

	virtual Vector 		gradient (const Vector& vec ) const{
		return VectorDifferential<VectorFunction>(*this).gradient(vec);
	}

	virtual Matrix<FT> 	hessian (const Vector& vec ) const {
		return Matrix<FT>();
	}
};

} // End of namespace COPT


#endif