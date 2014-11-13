//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef FUNC_REPOSITORY_H
#define FUNC_REPOSITORY_H

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
	typedef 	VectorFunction<VT>				Function;
	typedef 	typename Function::Vector 				Vector;
	typedef		typename Function::FT 					FT;
	Vector 			__weight;
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

template<class VT>
class TestQuadFunction
	: public VectorFunction<VT>
{
private:
	typedef VectorFunction<VT> Function;
	typedef 	typename Function::Vector 				Vector;
	typedef		typename Function::FT 					FT;

public:
	TestQuadFunction() {}
	FT operator() (const Vector& vec) const{
		return (vec[0]-50)*(vec[0]-50) + (vec[1]-10)*(vec[1]-20);
	}

};


template<class VT>
class TestQuadFunctionWithDiff
	: public VectorFunction<VT>
{
private:
	typedef VectorFunction<VT> Function;
	typedef 	typename Function::Vector 				Vector;
	typedef		typename Function::FT 					FT;

public:
	TestQuadFunctionWithDiff() {}
	FT operator() (const Vector& vec) const{
		return (vec[0]-50)*(vec[0]-50) + (vec[1]-10)*(vec[1]-20);
	}
	Vector gradient(const Vector& vec) const{
		Vector result(2);
		result[0] = 2*vec[0]-100;
		result[1] = 2*vec[1]-30;
		return result;
	}

};

template<class VT>
class RosenbrockFunction
	: public VectorFunction<VT>
{
public:
	RosenbrockFunction(){this->__dim = 2;}
	typename VectorFunction<VT>::ScalarType operator() (const VT& vec ) const{
		return (100*(vec[1]-vec[0]*vec[0])*(vec[1]-vec[0]*vec[0])+(1-vec[0])*(1-vec[0]));
	}
};

};

#endif