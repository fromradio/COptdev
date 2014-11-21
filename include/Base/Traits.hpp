// Copyright (C) 2014 Ruimin Wang ruimin.wang13@gmail.com
// Copyright (C) 2014 MathU

#ifndef COPT_TRAITS_H
#define COPT_TRAITS_H

namespace COPT
{
template<class T1>
struct get_pod_type
{ typedef T1 type;};

template<class T2>
struct get_pod_type<std::complex<T2> >
{ typedef T2 type;};

template<class T1>
struct is_float
{ static const bool value = false;};

template<>
struct is_float<float>
{ static const bool value = true;};

template<class T1>
struct is_double
{ static const bool value = false;};

template<>
struct is_double<double>
{ static const bool value = true;};

template<class T2>
struct is_real
{ static const bool value = false;};

template<>
struct is_real<float>
{ static const bool value = true;};

template<>
struct is_real<double>
{ static const bool value = true;};

template<class T1>
struct is_complex
{ static const bool value = false;};

template<class eT>
struct is_complex<std::complex<eT> >
{ static const bool value = false;};

template<>
struct is_complex<std::complex<float> >
{ static const bool value = true;};

template<>
struct is_complex<std::complex<double> >
{ static const bool value = true;};

template<class eT>
struct is_complex_float
{ static const bool value = false;};

template<>
struct is_complex_float<std::complex<float> >
{ static const bool value = true;};

template<class eT>
struct is_complex_double
{ static const bool value = false;};

template<>
struct is_complex_double<std::complex<double> >
{ static const bool value = true;};



template<class T>
struct is_scalar
{ static const bool value = 
					is_scalar<T>::value||
					is_complex<T>::value;
};



template<class T>
class VectorBase;
template<class T>
class MatrixBase;

/*		A trait class describing basic types that might be
 *		used in a numerical solver. A solver should take trait
 *		as template for flexibility.
 */
template<class T>
class KernelTrait
{
public:
	typedef T 					ScalarType;
	typedef VectorBase<T>		Vector;
	typedef MatrixBase<T>		Matrix;
};


/** traits of constraints and functions*/
struct linear_constraint_tag{};
struct quadratic_constraint_tag{};
struct non_linear_constraint_tag{};

template<class Constraint>
struct constraint_trait{
	typedef typename Constraint::constraint_category	constraint_category;
};

}// End of namespace COPT


#endif