// This file is part of COPT, a light-weight C++ based optimization open source library
//
// Copyright (C) 2015 Ruimin Wang <ruimin.wang13@gmail.com>
// Copyright (C) 2015 MathU
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#ifndef FUNC_REPOSITORY_HPP__
#define FUNC_REPOSITORY_HPP__

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
	VT gradient(const VT& vec) const{
		VT result(2);
		result[0] = 400*vec[0]*vec[0]*vec[0]-400*vec[0]*vec[1]+2*vec[0]-2;
		result[1] = 200*vec[1]-200*vec[0]*vec[0];
		return result;
	}
};

template<class T>
class VectorFunctionSystem
{
private:
	typedef 			VectorBase<T>			Vector;
	typedef 			MatrixBase<T>			Matrix;
	int __m;
	int __n;
public:
	VectorFunctionSystem()
		:
		__m(0),
		__n(0)
	{}
	VectorFunctionSystem(int m,int n)
		:
		__m(m),
		__n(n)
	{}
	~VectorFunctionSystem(){}

	Vector functionValue(const Vector& vec) const
	{
		Vector f(__m);
	//	1st
	//	f[0] = vec[0]*vec[0];
	//	f[1] = vec[1]*vec[1];
	//	2nd
	//	f[0] = (vec[0]-50)*(vec[0]-50) + (vec[1]-10)*(vec[1]-20);
	//	f[1] = (100*(vec[1]-vec[0]*vec[0])*(vec[1]-vec[0]*vec[0])+(1-vec[0])*(1-vec[0]));
	//	3rd
		f[0] = vec[0] + 10*vec[1]; 
		f[1] = sqrt(5)*(vec[2]-vec[3]);
		f[2] = (vec[1] - 2*vec[2]) * (vec[1] - 2*vec[2]);
		f[3] = sqrt(10) * (vec[0]-vec[3]) * (vec[0]-vec[3]);
		return f;
	}

	Matrix jacobiFunction(const Vector& vec) const
	{
		Matrix J(__m, __n);
	//	1st
	//	J(0, 0) = 2*vec[0];
	//	J(0, 1) = 0;
	//	J(1, 0) = 0;
	//	J(1, 1) = 2*vec[1];
	//	2nd
	//	J(0, 0) = 2*vec[0]-100;
	//	J(0, 1) = 2*vec[1]-30;
	//	J(1, 0) = 400*vec[0]*vec[0]*vec[0]-400*vec[0]*vec[1]+2*vec[0]-2;
	//	J(1, 1) = 200*vec[1]-200*vec[0]*vec[0];
	//	3rd
		J(0, 0) = 1;
		J(0, 1) = 10;
		J(1, 2) = sqrt(5);
		J(1, 3) = -sqrt(5);
		J(2, 1) = 2*(vec[1] - 2*vec[2]);
		J(2, 2) = -4*(vec[1] - 2*vec[2]);
		J(3, 0) = 2*sqrt(10)*(vec[0]-vec[3]);
		J(3, 3) = -2*sqrt(10)*(vec[0]-vec[3]);
		return J;

	}	

};

};

#endif
