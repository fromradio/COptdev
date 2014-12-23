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
struct is_size
{ static const bool value = false;};

// template<>
// struct is_size<size_t>
// { static const bool value = true;};

// template<>
// struct is_size<longsize>
// { static const boll value = true;};

// template<class T>
// struct integer_type
// {typedef T type;};
// template<>
// struct integer_type<size_t>
// { typedef int 	type;};

// template<>
// struct integer_type<longsize>
// { typedef COPTlong longsize;};

template<class T>
struct is_index
{ 
	typedef T size;
	static const bool value = false;
};

template<>
struct is_index<int>
{ 
	typedef unsigned int size;
	static const bool value = true;
};

template<>
struct is_index<long>
{ 
	typedef unsigned long size;
	static const bool value = true;
};

template<class T>
struct is_unsigned_size
{ 
	typedef T index;
	static const bool value = false;
};

template<>
struct is_unsigned_size<unsigned int>
{ 
	typedef 	int 			index;
	static const bool value = true;
};

template<>
struct is_unsigned_size<unsigned long>
{ 
	typedef 	long 			index;
	static const bool value = true;
};




template<class T>
struct is_scalar
{ static const bool value = 
					is_scalar<T>::value||
					is_complex<T>::value;
};


template<class T,class I>
class Array;
template<class T,class I>
class VectorBase;
template<class T,class I>
class MatrixBase;
template<class T,class I>
class SpMatrixBase;

/*		A trait class describing basic types that might be
 *		used in a numerical solver. A solver should take trait
 *		as template for flexibility.
 */
template<class T,class I = int >
class KernelTrait
{
public:
	typedef T 											scalar;
	typedef typename get_pod_type<scalar>::type 		podscalar;
	typedef I 											index;
	typedef typename is_index<I>::size					size;
	// whether the kernel is valid:
	static const bool valid  = is_scalar<T>::value&&is_index<I>::value;

	typedef Array<T,I> 									Array;
	typedef VectorBase<T,I>								Vector;
	typedef MatrixBase<T,I>								Matrix;
	typedef SpMatrixBase<T,I>							SpMatrix;
};

/** tags */
/** base tag for copt object */
struct copt_object{};
/** objects storing data */
struct data_object:copt_object{}; 				// the data object
struct referred_array:data_object{};			// refered array
struct array_object:data_object{}; 				// array
struct matrix_object:data_object{};				// normal matrix	
struct vector_object:data_object{};				// normal vector
struct sp_matrix_object:data_object{};			// sparse matrix

/** function objects */
struct scalar_func_object:copt_object{};		// scalar function
struct vector_func_object:copt_object{};		// vector function

/** problem object */
struct problem_object:copt_object{};			// problem
struct linear_programming:problem_object{};
struct lasso_problem:problem_object{};			// lasso problem

/** solver for problems */
struct solver_object:copt_object{};				// solver
struct linear_solver:solver_object{};			// linear solver
struct admm_solver:solver_object{};				// admm solver
struct proximal_solver:solver_object{}; 		// proximal solver
struct fista_solver:solver_object{};			// FISTA solver

/** constraint objects */
struct constraint_object:copt_object{};			// constraint

/** no time statistics tag */
struct time_stat_object:copt_object{};			// time statistics object
struct no_time_stat_tag:time_stat_object{};
struct solver_time_stat_tag:time_stat_object{};
/** traits of constraints and functions*/
struct linear_constraint:constraint_object{};
struct linear_leq_constraint:linear_constraint{};
struct linear_neq_constraint:linear_constraint{};
struct linear_eq_constraint:linear_constraint{};
struct quadratic_constraint_tag:constraint_object{};
struct non_linear_constraint_tag:constraint_object{};

template<class Constraint>
struct constraint_trait{
	typedef typename Constraint::constraint_category	constraint_category;
};



/** trais of functions */
struct scalar_function_tag{};
struct log_scalar_function_tag:public scalar_function_tag{};
struct abs_scalar_function_tag:public scalar_function_tag{};


/** force assignment from complex to real */
template<class T>
void ForceAssignment(const T& t1,std::complex<T>& t2)
{
	// real to complex
	if( is_real<T>::value )
		t2 = std::complex<T>(t1,0);
}

template<class T>
void ForceAssignment(const std::complex<T>&t1 , T& t2 )
{
	if ( is_real<T>::value )
		t2 = t1.real();
}

template<class T>
void ForceAssignment(const T& t1, T& t2)
{
	t2 = t1;
}

}// End of namespace COPT


#endif