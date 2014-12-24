//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef LP_PROBLEM_HPP__
#define LP_PROBLEM_HPP__


namespace COPT
{


/*		A base design for general constraint. An optimization constraint 
 *		constains several parts like leq constraint, equal constraint etc.
 */

template<class kernel>
class GeneralConstraint
{

	typedef typename kernel::Vector 		Vector;
public:
	/*		The type of constraints. Typically the type can be reduced into
	 *		three: equal, strictly less than and less than.
	 */
	enum ConType{
		LeqStrict,
		Leq,
		Eq,
		Neq,
		NeqStrict,
		LowerBound,
		UpperBound,
		BothBound
	};

private:
	/** the type of the constraint */
	ConType 			__type;

public:
	typedef constraint_object 				ObjectCategory;

	GeneralConstraint( const char c , const char s );
	virtual ~GeneralConstraint(){}
	virtual bool feasible( const Vector& x ) const = 0;

	/** getter and setter */
	//%{
	/** set the constraint type */
	void setConstraintType(const char c , const char s = 'n' );
	ConType constraintType() const;
	//%}
};

/****************Implementation of class 'GeneralConstraint'***************/
template<class kernel>
GeneralConstraint<kernel>::GeneralConstraint( const char c , const char s )
{
	setConstraintType(c,s);
}

template<class kernel>
void GeneralConstraint<kernel>::setConstraintType(const char c , const char s)
{
	if ( c == 'l' || c == 'L' )
	{
		if ( s == 'n' || s == 'N' )
			__type = Leq;
		else if ( s == 'y' || s == 'Y' )
			__type = LeqStrict;
		else
		{
			std::cerr<<"Warning: the input for strict or not should be 'n' or 'y'. 'y' would be taken for other input!"<<std::endl;
			__type = LeqStrict;
		}
	}
	else if ( c == 'n' || c == 'N' )
	{
		if ( s == 'n' || s =='N' )
			__type = Neq;
		else if ( s == 'y' || s == 'Y' )
			__type = NeqStrict;
		else
		{
			std::cerr<<"Warning: the input for strict or not should be 'n' or 'y'. 'y' would be taken for other input!"<<std::endl;
			__type = NeqStrict;
		}
	}
	else if ( c == 'e' || c == 'E' )
	{
		__type = Eq;
	}
	else if ( c == 'b' || c == 'B' )
	{
		if ( s == 'l' || s == 'L' )
			__type = LowerBound;
		else if ( s == 'u' || s == 'U' )
			__type = UpperBound;
		else if ( s == 'b' || s == 'B' )
			__type = BothBound;
	}
	else
	{
		throw COException("Wrong type for type of linear constraint. The input for the constraint should be 'l' or 'L' for leq, 'n' or 'N' for neq and 'e'or 'E' for eq!");
	}
}

template<class kernel>
typename GeneralConstraint<kernel>::ConType GeneralConstraint<kernel>::constraintType() const
{
	return __type;
}

template<class kernel>
class BoundConstraint
	:
	public GeneralConstraint<kernel>,
	noncopyable
{
private:
	typedef GeneralConstraint<kernel> 		Base;
	typedef typename kernel::scalar 		scalar;
	typedef typename kernel::index 			index;
	typedef typename kernel::Vector 		Vector;
	typedef typename kernel::Matrix 		Matrix;

	/** the lower bound */
	scalar 			__l;
	/** the upper bound */
	scalar 			__u;

public:

	/** constructor and deconstructor */
	//%{
	/*		Default constructor assumes that the boundary constraint
	 *		is both side with 0<=x<=0
	 */
	BoundConstraint();
	BoundConstraint(
		const scalar sc,
		const char s = 'l');
	BoundConstraint(
		const scalar l,
		const scalar u);
	//%}

	/** check the feasibility of the constraint */
	bool feasible ( const Vector& x ) const;

	/** getter and setter */
	//%{
	void setLowerBound( const scalar l );
	scalar lowerBound() const;
	void setUpperBound( const scalar u );
	scalar upperBound() const;
	//%}
};

/***************Implementation of class 'BoundConstraint'************/
template<class kernel>
BoundConstraint<kernel>::BoundConstraint()
	:
	Base('b','b'),
	__l(0.0),
	__u(0.0)
{
}

template<class kernel>
BoundConstraint<kernel>::BoundConstraint(
	const scalar sc,
	const char s)
	:
	Base('b',s)
{
	if ( s == 'l' )
		__l = sc;
	else if ( s == 'u' )
		__u = sc;
	else
		throw COException("Unknown type for boundary constraint!");
}

template<class kernel>
BoundConstraint<kernel>::BoundConstraint(
	const scalar l,
	const scalar u)
	:
	Base('b','b'),
	__l(l),
	__u(u)
{
}

template<class kernel>
bool BoundConstraint<kernel>::feasible( const Vector& x ) const
{
	switch (this->constraintType())
	{
		case Base::BothBound:
		{
			return (x>=__l&&x<=__u);
		}
		break;
		case Base::LowerBound:
		{
			return (x>=__l);
		}
		break;
		case Base::UpperBound:
		{
			return (x<=__u);
		}
		break;
		default:
		{
			throw COException("Error: Unknown boundary constraint type!");
		}
		break;   
	}
}

template<class kernel>
void BoundConstraint<kernel>::setLowerBound( const scalar l )
{
	__l = l;
}

template<class kernel>
typename BoundConstraint<kernel>::scalar BoundConstraint<kernel>::lowerBound() const
{
	return __l;
}

template<class kernel>
void BoundConstraint<kernel>::setUpperBound( const scalar u )
{
	__u = u;
}

template<class kernel>
typename BoundConstraint<kernel>::scalar BoundConstraint<kernel>::upperBound() const
{
	return __u;
}

//////////////End of implementation of 'BoundConstraint'


/*		A description of a linear constraint. 
 *
 *
 */
template<class kernel>
class LinearConstraint
	:
	public GeneralConstraint<kernel>,
	noncopyable
{
private:
	typedef GeneralConstraint<kernel> 		Base;
	typedef typename kernel::scalar 		scalar;
	typedef typename kernel::index 			index;
	typedef typename kernel::Vector 		Vector;
	typedef typename kernel::Matrix 		Matrix;

public:

protected:

	/* 	whether the constraint contains a matrix 
	 *	if __is_matrix is true, the left hand of constraint has form Ax.
	 * 	if __is_matrix is false, the left hand of constraint is x.
	 */
	bool 				__is_matrix;

	/* 	the left hand matrix */
	Matrix 				__A;
	/* 	the right hand vector */
	Vector 				__b;
	/* 	the dimension of the constraint, i.e. the dimension of x. */
	index 				__dim;

public:
	typedef linear_constraint 				ObjectCategory;

	/** deconstructor and constructor */
	//%{
	LinearConstraint(
		const index dim,
		const char c = 'l',
		const char s = 'n'
		);
	LinearConstraint(
		const Vector& b,
		const char c = 'l',
		const char s = 'n'
		);
	LinearConstraint(
		const Matrix& A,
		const Vector& b,
		const char c = 'l',
		const char s = 'n'
		);
	virtual ~LinearConstraint(){}
	//%}

	

	/** whether the input is a feasible point */
	bool feasible(const Vector& x) const;

	/** getter and setter */
	//%{
	void setMatA( const Matrix& A );
	const Matrix& matA() const;
	void setRhB( const Vector& b );
	const Vector& rhB() const;
	//%}
};

/*************Implementation of class 'LinearConstraint'**************/
template<class kernel>
LinearConstraint<kernel>::LinearConstraint(
	const index dim,
	const char c,
	const char s)
	:
	GeneralConstraint<kernel>(c,s),
	__is_matrix(false),
	__b(dim),
	__dim(dim)
{
}

template<class kernel>
LinearConstraint<kernel>::LinearConstraint(
	const Vector& b,
	const char c,
	const char s)
	:
	GeneralConstraint<kernel>(c,s),
	__is_matrix(false),
	__b(b),
	__dim(b.size())
{
}

template<class kernel>
LinearConstraint<kernel>::LinearConstraint(
	const Matrix& A,
	const Vector& b,
	const char c,
	const char s)
	:
	GeneralConstraint<kernel>(c,s),
	__is_matrix(true),
	__A(A),
	__b(b),
	__dim(A.cols())
{
}

template<class kernel>
void LinearConstraint<kernel>::setMatA(const Matrix& A)
{
	__A = A;
}

template<class kernel>
const typename LinearConstraint<kernel>::Matrix& LinearConstraint<kernel>::matA() const
{
	return __A;
}

template<class kernel>
void LinearConstraint<kernel>::setRhB(const Vector& b)
{
	__b = b;
}

template<class kernel>
const typename LinearConstraint<kernel>::Vector& LinearConstraint<kernel>::rhB() const
{
	return __b;
}

template<class kernel>
bool LinearConstraint<kernel>::feasible( const Vector& x ) const
{
	switch (this->constraintType())
	{
	case Base::Leq:
	{
		if ( __is_matrix )
			return (__A*x<=__b);
		else
			return x<=__b;
	}
	break;
	case Base::LeqStrict:
	{
		if ( __is_matrix )
			return (__A*x<__b);
		else
			return x<__b;
	}
	break;
	case Base::Eq:
	{
		if ( __is_matrix )
			return (__A*x==__b);
		else
			return x==__b;
	}
	break;
	case Base::NeqStrict:
	{
		if ( __is_matrix )
			return (__A*x>__b);
		else
			return x>__b;
	}
	break;
	case Base::Neq:
	{
		if ( __is_matrix )
			return (__A*x>=__b);
		else
			return x>=__b;
	}
	break;
	default:
	{
		throw COException("Unknown type for linear constraint!");
	}
	break;
	}
}

///////////////////End of implementation of 'LinearConstraint'


/*		The description for linear programming problems. 
 *		An LP problem contains a linear objective and several linear constraints.
 *		A linear objective function of variable 'x' is written as w^Tx where
 *		w is a coefficient vector. Linear constraints contain several part as
 *		Ax<=b, Bx=0, boundary, x>=0 and so on.
 *		A standard form of LP problem is:
 *						minimize w^Tx
 *							subject to Ax<=b, x>=0
 *		Note that any other LP problem can be transformed into a standard form.
 *		Typical solvers also assume that LP problem is written in a standrard form.
 * 		In standard linear programming problem, the formulation is very simple.
 *		One is the weight vector w and one is the Ax<=b constraint. The constraint
 * 		that x>=0 is not recorded but taken as default.
 */

template<class kernel>
class LPProblem
	:
	public VectorProblem<kernel>
{

	typedef typename kernel::scalar 		scalar;
	typedef typename kernel::index 			index;
	typedef typename kernel::Vector 		Vector;
	typedef typename kernel::Matrix 		Matrix;

	/** the weight vector */
	Vector 						__w;
	/** the Ax<=b constraint*/
	LinearConstraint<kernel> 	__axb;

public:

	/** constructor and deconstructor */
	//%{
	LPProblem( const Vector& w , const Matrix& A , const Matrix &b );
	~LPProblem(){}
	//%}

	/** compute objective function of LP problem */
	scalar objective( const Vector& x ) const;

	/** valid */
	bool isValidInput( const Vector& x ) const;

	/** validation of the problem */
	bool isValid( ) const;
	
	/** whether one input is feasible */
	bool feasible ( const Vector& x ) const;

	

};

/**************Implementation of class 'LPProblem'**************/
template<class kernel>
LPProblem<kernel>::LPProblem( const Vector& w , const Matrix& A , const Matrix &b )
	:
	VectorProblem<kernel>(A.cols()),
	__w(w),
	__axb(A,b)
{
	if (!isValid())
	{
		std::cerr<<"Linear programming problem error: the problem is not valid! You should better check it more carefully!"<<std::endl;
	}
}

template<class kernel>
typename LPProblem<kernel>::scalar LPProblem<kernel>::objective(const Vector& x )const
{
	if ( x.size() != __w.size() )
		throw COException("Linear programming objective function computation error: the size of weight function and x are not consistent!");
	return __w.dot(x);
}

template<class kernel>
bool LPProblem<kernel>::isValidInput( const Vector& x ) const
{
	return x.size()==this->dimension();
}

template<class kernel>
bool LPProblem<kernel>::isValid( ) const
{
	if (__w.size() != __axb.matA().cols() )
		return false;
	else
		return __axb.matA().rows()==__axb.rhB().size();
}

template<class kernel>
bool LPProblem<kernel>::feasible(const Vector& x ) const
{
	/** Ax<=b && x>=0 */
	return (__axb.feasible(x)&&(x>=0));
}

//////////////////////End of implementation of class 'LPProblem'


}
#endif