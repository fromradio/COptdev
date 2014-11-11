// Function.h
#ifndef FUNCTION_H
#define FUNCTION_H

#include "Vector.h"
#include "Matrix.h"
// #include "BasicMath.h"
// #include "FuncPara.h"
#include "NumDifferential.h"

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
	typedef 		T 			 FT;
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
	typedef 		VT 						Vector;
	typedef 		typename Vector::FT 	FT;

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

};

/*
	parameters of functions
		a double vector with a few functionalities
*/


/*
	a general parent class for function
*/
// a general class for scalar function
// class ScalarFunc{
// public:
// 	// function evaluation
// 	virtual 	double operator()	( double x ) = 0;
// 	// function differential
// 	virtual 	double diff 		( double x ) = 0;
// 	ScalarFunc(){}
// protected:
// };

// class Function
// {
// public:
// 	// Virtual operations for functions
// 	/*
// 		evaluation for a vector
// 	*/
// 	virtual double operator()(const Vector& vec) = 0;
// 	/*
// 		evaluation for differential if possible
// 	*/
// 	bool diffable() {return __exact_diff;}
// 	virtual void diff(Vector& value,const Vector& vec) = 0;

// 	Function();
// protected:
// 	bool	__exact_diff;
// };
// /*
// 	Vector Function
// */
// class VectorFunc
// {
// public:
// 	virtual Vector 	operator() 	( const Vector& vec )					= 0;
// 	virtual void 	evaluation 	( Vector& value , const Vector& vec )	= 0;
// 	virtual void 	diff 		( Matrix& value , const Vector& vec )	= 0;
// protected:
// };

// /*
// 	Differential function of 
// */


// /*
// 	an example for 'CosineFunction'
// */
// class ScalarCosineFunc:public ScalarFunc
// {
// public:
// 	double 		operator()	( double x );
// 	double 		diff 		( double x );
// 	ScalarCosineFunc		( double l = 1.0 );
// private:
// 	FuncPara<double> __fp;
// };
// /*
// 	an example for 'Cosine Function'
// */
// class CosineFunc:public Function
// {
// public:
// 	double 		operator()	(const Vector& vec);
// 	void 		diff 		(Vector& value,const Vector& vec);
// 	CosineFunc 				(const FuncPara<double>& fp);
// };
// /*
// 	an example for 'Vector Cosine Function'
// 	y_i = cos(a_i x_i)
// */
// class VectorCosineFunc : public VectorFunc
// {
// public:
// 	Vector 			operator()	( const Vector& vec );
// 	void 			evaluation	( Vector& value , const Vector& vec ) {}
// 	void			diff 		( Matrix& value , const Vector& vec ) {}
// 	VectorCosineFunc			( ) {}
// 	VectorCosineFunc			( int dim , double* vec ){}
// private:
// 	int							__dim;
// };


#endif