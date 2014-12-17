// 		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef BP_PROBLEM_HPP__
#define BP_PROBLEM_HPP__

namespace COPT
{

/*		Basis pursuit problem is a typical and traditional problem
 *		in optimization field. The model is that:
 *			minimize \|x\|_1
 *				subject to Ax = b
 */
template<class kernel>
class BPProblem
{
private:
	typedef typename kernel::scalar 				scalar;
	typedef typename kernel::index 					index;
	typedef typename kernel::Vector 				Vector;
	typedef typename kernel::Matrix 				Matrix;


	LinearConstraint<kernel> 		__axb;


public:

	BPProblem( const Matrix& A , const Vector& b );
	scalar objective ( const Vector& x ) const;

	/** judge whether the input x is feasible or not */
	bool feasible( const Vector& x ) const;

	/** judge whether the input x is valid or not */
	bool isValidInput(const Vector& x ) const;

	/** judge whether the problem is valid */
	bool isValid() const;

};

/**************Implementation of class 'BPProblem'*****************/
template<class kernel>
BPProblem<kernel>::BPProblem( const Matrix& A , const Vector& b )
	:
	VectorProblem<kernel>(A.cols()),
	__axb(A,b)
{
}

template<class kernel>
bool BPProblem<kernel>::isValidInput( const Vector& x ) const
{
	if ( x.size() == __axb.matA().cols() )
		return true;
	else
		return false;
}

template<class kernel>
typename BPProblem<kernel>::scalar BPProblem<kernel>::objective(const Vector& x) const
{
	return x.absNorm();
}

template<class kernel>
bool BPProblem<kernel>::feasible( const Vector& x ) const
{
	return __axb.feasible(x);
}

template<class kernel>
bool BPProblem<kernel>::isValid( ) const
{
	return __axb.matA().rows()==__axb.rhB().size();
}

///////////////End of implementation of 'BPProblem'


}// End of namespace COPT

#endif