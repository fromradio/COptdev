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

template<class T,class scalar>
T computeProximal( const LogFunction& func , const T& v , const scalar lambda , const log_scalar_function_tag& )
{
	return (v+std::sqrt(v*v+4*lambda))/2;
}

template<class T,class scalar>
T computeProximal( const AbsFunction& func , const T& v , const scalar lambda , const abs_scalar_function_tag& )
{
	if ( v >= lambda )
		return (v-lambda);
	else if ( v <= -lambda )
		return (v+lambda);
	else
		return 0;
}

template<class Function,class T,class scalar>
T computeProximal( const Function& func , const T& v , const scalar lambda )
{
	return computeProximal(func,v,lambda,typename Function::function_category());
}
}

#endif