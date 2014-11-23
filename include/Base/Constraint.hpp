// Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
// Copyright (C) MathU

namespace COPT
{

/*		Basic class for describing general constraint. The class
 *		takes the type of data as template. A virtual function judging
 *		whether the input data is feasible is defined. Every constraint
 *		should derive the method
 */	
template<class kernel>
class Constraint
{
	typedef typename	kernel::Vector 			Vector;
public:
	//%{
	/** Default constructor*/
	Constraint(){}
	/** Deconstructor*/
	virtual ~Constraint(){}
	//%}

	/** Judge whether the input is feasible*/
	virtual bool feasible(const Vector& ) const = 0;
};


template<class kernel>
class LinearEqualConstraint
	:
	Constraint<VectorBase<kernel> >
{
private:
	typedef typename	kernel::Vector 			Vector;
	typedef typename 	kernel::Matrix 			Matrix;
	typedef linear_constraint_tag 				constraint_category;
	/** variables*/
	//%{
	/** Matrix A*/
	Matrix 			__A;
	/** Vector b*/
	Vector 			__b;
	//%}
public:
	/** Constructor and Deconstructor*/
	//%{
	LinearEqualConstraint(const Matrix& A,const Vector& b);
	//%}

	bool feasible(const Vector& ) const;

	/** getter */
	//%{
	/** get the A matrix */
	const Matrix& mat() const;
	/** get the b vector */
	const Vector& vec() const;
	//%}
};

// Implementation
template<class kernel>
LinearEqualConstraint<kernel>::LinearEqualConstraint( const typename kernel::Matrix&A,const typename kernel::Vector&b)
	:
	__A(A),
	__b(b)
{
}

template<class kernel>
bool LinearEqualConstraint<kernel>::feasible( const typename kernel::Vector& x)const
{
	return (A*x==b);
}

template<class kernel>
const typename kernel::Matrix& LinearEqualConstraint<kernel>::mat() const
{
	return __A;
}

template<class kernel>
const typename kernel::Vector& LinearEqualConstraint<kernel>::vec() const
{
	return __b;
}
}