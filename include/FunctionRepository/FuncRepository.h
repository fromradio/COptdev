/*
	Powered by ‘MathU'
		Copyright@MathU
*/

#ifndef FUNC_REPOSITORY_H
#define FUNC_REPOSITORY_H

#include "Function.h"

namespace COPT
{

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


/*
 *					Pre-defined vector function examples
 *
 */
template<class VT>
class VectorCosineFunction
	: public VectorFunction<VT>
{
private:
	typedef 	typename VectorFunction<VT>::Vector 	Vector;
	typedef		typename VectorFunction<VT>::FT 		FT;
	VT 			__weight;
public:
	VectorCosineFunction(const Vector& w):__weight(w){this->__dim = w.size();}
	FT operator() (const Vector& vec) const {
		if (vec.size() != this->__dim ){
			// std::cout<<vec.size()<<' '<<this->__dim<<std::endl;
			throw COException("Function error: dimension must be the same!");
		}
		return cos(__weight.dot(vec));
	}
};

};

#endif