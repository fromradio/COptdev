// Function.h
#ifndef FUNCTION_H
#define FUNCTION_H

#include "Vector.h"
#include "Matrix.h"
// #include "BasicMath.h"
// #include "FuncPara.h"

/*
	Second version of 'Functions' classes
*/
namespace COPT{
/*
	Scalar function
		input: scalar x
		output: scalar y
	template
		T stands for the float type that is used, for example double or float

*/
template<class T>
class ScalarFunction{
protected:
	// whether the function has true differential
	bool		__truediff;
public:
	typedef T FT;
	ScalarFunction (  ):__truediff(false){}
	virtual ~ScalarFunction() {}
	// basic operation computing the value of the function
	virtual FT operator() 	( FT x ) const = 0;
	// basic operation computing the differential of the function
	virtual FT diff  		( FT x ) const = 0;
	// whethre the function has accurate difference
	bool		diffable	() const {return __truediff;}
};

template <class FT>
class CosineFunction:public ScalarFunction<FT>{
private:
	FT __l;
public:
	CosineFunction( FT l ):__l(l){}
	FT operator() ( FT x ) const {return cos(__l*x);}
	FT diff(FT x) const {return -__l*sin(__l*x);}
};

template<class FT>
class LinearFunction:public ScalarFunction<FT>{
private:
	FT __a;
	FT __b;
public:
	LinearFunction(FT a,FT b):__a(a),__b(b){}
	FT operator() ( FT x ) const {return __a*x+__b;}
	FT diff(FT x) const {return __a;}
};

template<class FT>
class QuadFunction:public ScalarFunction<FT>{
private:
	FT __a;
	FT __b;
public:
	QuadFunction(FT a,FT b):__a(a),__b(b){}
	FT operator() ( FT x ) const {return cos(x)-x*x*x;}
	FT diff(FT x) const {return -sin(x)-3*x*x;}
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