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
	:
	public VectorProblem<kernel>
{
private:
	typedef typename kernel::scalar 				scalar;
	typedef typename kernel::index 					index;
	typedef typename kernel::Vector 				Vector;
	typedef typename kernel::Matrix 				Matrix;


	LinearConstraint<kernel> 		__axb;


public:

	/** constructor and deconstructor */
	//%{
	BPProblem(
		const index dim
		);
	BPProblem( 
		const index dim, 		// the dimension of the problem
		const Matrix& A, 		// the A matrix in Ax=b
		const Vector& b 		// the right hand vector b
		);
	//%}

	scalar objective ( const Vector& x ) const;

	/** judge whether the input x is feasible or not */
	bool feasible( const Vector& x ) const;

	/** judge whether the input x is valid or not */
	bool isValidInput(const Vector& x ) const;

	/** judge whether the problem is valid */
	bool isValid() const;

	/** getter and setter */
	//%{
	void setMatA( const Matrix& A );
	const Matrix& matA( ) const;
	void setRhB( const Vector& b);
	const Vector& rhB( ) const;
	//%}

};

/**************Implementation of class 'BPProblem'*****************/

template<class kernel>
BPProblem<kernel>::BPProblem(
	const index dim
	)
	:
	VectorProblem<kernel>(dim)
{
}

template<class kernel>
BPProblem<kernel>::BPProblem( 
	const index dim,
	const Matrix& A , 
	const Vector& b )
	:
	VectorProblem<kernel>(dim),
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
	if ( this->dimension() != __axb.matA().cols() )
		return false;
	return __axb.matA().rows()==__axb.rhB().size();
}

template<class kernel>
void BPProblem<kernel>::setMatA( const Matrix& A )
{
	__axb.setMatA(A);
}

template<class kernel>
const typename BPProblem<kernel>::Matrix& BPProblem<kernel>::matA() const
{
	return __axb.matA();
}

template<class kernel>
void BPProblem<kernel>::setRhB(const Vector& b)
{
	__axb.setRhB(b);
}

template<class kernel>
const typename BPProblem<kernel>::Vector& BPProblem<kernel>::rhB() const
{
	return __axb.rhB();
}

///////////////End of implementation of 'BPProblem'

template<class Problem,class Time = NoTimeStatistics>
class BPSolver
	:
	public GeneralSolver<typename Problem::KernelTrait,Time>,
	noncopyable
{
private:

	typedef typename Problem::KernelTrait				kernel;
	typedef typename kernel::scalar 					scalar;
	typedef typename kernel::index 						index;
	typedef typename kernel::Vector 					Vector;
	typedef typename kernel::Matrix 					Matrix;

	/** the problem */
	const Problem& 					__p;
	/** storing the result */
	Vector 							__x;
	/** y vector */
	Vector 							__y;
	/** storing the matrix for y sub problem */
	Matrix 							__mat;
	/** mu */
	scalar 							__mu;
	/** lambda for x */
	Vector 							__lambda_x;
	/** lambda for y */
	Vector 							__lambda_y;
	/** the linear solver */
#ifdef EIGEN
	Eigen::LDLT<Eigen::MatrixXd>	__ldlt;
#endif

	/** the type of linear solver that is used */
	LinearSolverType 				__linear_type;

	/** the pointer to the linear solver */
	LinearSolver<Matrix> 			*__linear_solver;

	/** using the LU factorization for storing the nessessary information */
	LU<Matrix>						__lu;

	/** before solving */
	void solvingBegin();

	/*		compute the problem, the decomposition of A^TA+I is computed at first
	 *
	 */
	void doCompute();

	scalar doOneIteration();

	/** solve x sub problem */
	void xSubproblem();

	/** solve y sub problem */
	void ySubproblem();

	/** update of lambda x */
	void lambdaXUpdate();

	/** update of lambda y */
	void lambdaYUpdate();

	/** generate solver */
	void genreateLinearSolver();

public:

	/** constructor and deconstructor */
	//%{
	BPSolver( 
		const Problem& p,
		const index maxiteration = 10000,
		const scalar thresh = 1e-10 ,
		const enum LinearSolverType type = COPT::CholeskySolver );
	~BPSolver();
	//%}

	/** setter and getter */
	//%{
	/** compute the objective value */
	scalar objective() const;
	/** obtain the result */
	const Vector& result() const;
	/** set the linear solve type */
	void setLinearSolverType(const enum LinearSolverType type);
	//%}


};

/*************Implementation of class 'BPSolver'******************/

template<class Problem,class Time>
void BPSolver<Problem,Time>::solvingBegin()
{
	__x.resize(__p.dimension());
	__y.resize(__p.dimension());
	__lambda_x.resize(__p.dimension());
	__lambda_y.resize(__p.rhB().size());
}

template<class Problem,class Time>
void BPSolver<Problem,Time>::doCompute()
{
	// get the matrix A^TA+I
	Matrix mtm;
	__p.matA().mtm(mtm);
	index dim = __p.dimension();
	for ( int i = 0 ; i < dim ; ++ i )
		mtm(i,i) += 1.0;
	SAFE_DELETE(__linear_solver);
	genreateLinearSolver();
	__linear_solver->compute(mtm);
	__mat = __p.matA().transpose();
// #ifdef EIGEN
// 	Eigen::MatrixXd m(dim,dim);
// 	for ( int i = 0 ; i < dim ; ++ i )
// 		for ( int j = 0 ; j < dim ; ++ j )
// 			m(i,j) = __mat(i,j);
// 	__ldlt.compute(m);
// #endif
}

template<class Problem,class Time>
typename BPSolver<Problem,Time>::scalar BPSolver<Problem,Time>::doOneIteration()
{
	xSubproblem();
	ySubproblem();
	lambdaXUpdate();
	lambdaYUpdate();

	return (__x-__y).squaredNorm();
}

template<class Problem,class Time>
void BPSolver<Problem,Time>::xSubproblem()
{
	computeProximal(AbsFunction(),__y+__mu*__lambda_x,__mu,__x);
}

template<class Problem,class Time>
void BPSolver<Problem,Time>::ySubproblem()
{
	__linear_solver->solve(__mat*(__p.rhB()-__mu*__lambda_y));
	// __y = __lu.solve(__mat*(__p.rhB()-__mu*__lambda_y));
// #ifdef EIGEN
// 	Vector rhb = __x-__mu*__lambda_x+__p.matA().transpose()*(__p.rhB()-__mu*__lambda_y);
// 	Eigen::VectorXd r(rhb.size());
// 	int dim = __p.dimension();
// 	for ( int i = 0 ; i < dim ; ++ i )
// 		r(i) = rhb(i);
// 	Eigen::VectorXd y = __ldlt.solve(r);
// 	for ( int i = 0 ; i< dim ; ++ i )
// 		__y(i) = y(i);
// #endif
}

template<class Problem,class Time>
void BPSolver<Problem,Time>::lambdaXUpdate()
{
	__lambda_x = __lambda_x + (1.0/__mu)*(__y-__x);
}

template<class Problem,class Time>
void BPSolver<Problem,Time>::lambdaYUpdate()
{
	__lambda_y = __lambda_y + (1.0/__mu)*(__p.matA()*__y-__p.rhB());
}

template<class Problem,class Time>
void BPSolver<Problem,Time>::genreateLinearSolver()
{
	SAFE_DELETE(__linear_solver);
	switch(__linear_type)
	{
		case LUSolver:
		{
			__linear_solver = new LU<Matrix>;
		}
		break;
		case QRSolver:
		{
			__linear_solver = new QR<Matrix>;
		}
		break;
		case CholeskySolver:
		{
			__linear_solver = new Cholesky<Matrix>;
		}
		break;
		case EigenWrap:
		{
			__linear_solver = new EigenSolver<Matrix>;
		}
		break;
		default:
		{
		throw COException("Unknown linear solver for BP solver!");
		}
		break;
	}
}

template<class Problem,class Time>
BPSolver<Problem,Time>::BPSolver( 
	const Problem& p,
	const index maxiteration,
	const scalar thresh,
	const enum LinearSolverType type
	 )
	:
	GeneralSolver<kernel,Time>(maxiteration,thresh),
	__p(p),
	__mu(1.0),
	__linear_type(type),
	__linear_solver(NULL)
{
	this->compute();
}

template<class Problem,class Time>
BPSolver<Problem,Time>::~BPSolver()
{
	SAFE_DELETE(__linear_solver);
}

template<class Problem,class Time>
typename BPSolver<Problem,Time>::scalar BPSolver<Problem,Time>::objective() const
{
	return __p.objective(__x);
}

template<class Problem,class Time>
const typename BPSolver<Problem,Time>::Vector& BPSolver<Problem,Time>::result() const
{
	return __x;
}

template<class Problem,class Time>
void BPSolver<Problem,Time>::setLinearSolverType( const enum LinearSolverType type )
{
	__linear_type = type;
}


}// End of namespace COPT

#endif