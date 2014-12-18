// 		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef PROXIMAL_OPERATORS_HPP__
#define PROXIMAL_OPERATORS_HPP__

/** 	Fast computation of famous proximal operators 
 */

namespace COPT
{

/*		calculation of proximal operators
 */
template<class Function>
class ProximalOperator
{

};

class QuardraticProximal
{
private:
public:

};

class LogFunction
{
public:
	typedef		log_scalar_function_tag 	function_category;
	template<class FT>
	FT operator()(const FT x){
		return -log(x);
	}
};

class AbsFunction
{
public:
	typedef 	abs_scalar_function_tag 	function_category;
	template<class FT>
	FT operator()(const FT x){
		return std::abs(x);
	}
};

template<class scalar>
scalar computeProximal( const LogFunction& func , const scalar v , const scalar lambda , const log_scalar_function_tag& )
{
	return (v+std::sqrt(v*v+4*lambda))/2;
}

template<class scalar>
scalar computeProximal( const AbsFunction& func , const scalar v , const scalar lambda , const abs_scalar_function_tag& )
{
	if ( v >= lambda )
		return (v-lambda);
	else if ( v <= -lambda )
		return (v+lambda);
	else
		return 0;
}

template<class Vector,class scalar>
void computeProximal( const AbsFunction& func, const Vector& v , const scalar lambda , const abs_scalar_function_tag& , Vector& x )
{
	if (x.size() != v.size() )
		x.resize(v.size());
	for ( int i = 0 ; i < v.size() ; ++ i )
	{
		x(i) = computeProximal(func,v(i),lambda);
	}
}

template<class Function,class scalar>
scalar computeProximal( const Function& func , const scalar v , const scalar lambda )
{
	return computeProximal(func,v,lambda,typename Function::function_category());
}

template<class Function,class Vector,class scalar>
void computeProximal( const Function& func , const Vector& v , const scalar lambda , Vector& x )
{
	computeProximal(func,v,lambda,typename Function::function_category(),x);
}
}

#endif