//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef LASSO_HPP__
#define LASSO_HPP__

namespace COPT
{

/*		A general design for solver. The solver derives from a Time
 *		Stastistics class to help it compute time cost of the solver.
 *
 */	
template<class Time=NoTimeStatistics>
class GeneralSolver
	:
	public Time
{
protected:
	virtual void doSolve() = 0;
	virtual void doCompute() = 0;
public:

	/** virtual deconstructor */
	virtual ~GeneralSolver(){}
	/** kernel function: solve the problem */
	void solve();
	/** kernel function: compute or initialize the problem */
	void compute();
};

template<class Time>
void GeneralSolver<Time>::compute()
{
	this->computationBegin();
	this->doCompute();
	this->computationEnd();
}

template<class Time>
void GeneralSolver<Time>::solve()
{
	this->solvingBegin();
	this->doSolve();
	this->solvingEnd();
}

/*		Proximal solver for lasso problem
 *
 *
 */
template<class Problem,class Time = NoTimeStatistics>
class LassoProximalSolver
	:
	public GeneralSolver<Time>,
	noncopyable
{
private:
	typedef typename Problem::KernelTrait::index		index;
	typedef typename Problem::KernelTrait::scalar 		scalar;
	typedef typename Problem::KernelTrait::Vector 		Vector;
	typedef typename Problem::KernelTrait::Matrix 		Matrix;

	/** private variables */
	//%{
	/** the reference to the problem */
	const Problem& 								__s;
	/** A^TA matrix */
	Matrix 										__mtm;
	/** A^T */
	Matrix 										__mt;
	/** the proximal parameter */
	scalar 										__mu;
	/** scaling factor of proximal parameter */
	scalar 										__beta;
	/** maximal iteration number */
	index 										__maxiteration;
	/** the result */
	Vector 										__x;

	Eigen::LDLT<Eigen::MatrixXd>				__linear_solver;
	//%}

	// void init();

	LassoProximalSolver();

	void doCompute();
	void doSolve();
public:

	LassoProximalSolver(
		const Problem& s , 
		const scalar mu = 0.5 , 
		const scalar beta = 1.0,
		const index maxiteration = 10000);

	/** compute the gradient of f */
	Vector fGradient(const Vector&);

	/** find the proximal of g */
	Vector gProximal(const Vector &);

	/** the f function */
	scalar f(const Vector& );

	/** the f_mu function */
	scalar fMu(const Vector&,const Vector&);

	/** a single iteraton of the solver */
	bool iteration(const Vector& xk,Vector& xn);

	/** setter and getter */
	//%{
	/** the final result of proximal solver */
	Vector result() const;
	//%}

};

/*		ADMM solver for solving lasso problem
 *
 *
 */
template<class Problem,class Time = NoTimeStatistics>
class LassoADMMSolver
	:
	public GeneralSolver<Time>
{
private:
	typedef typename Problem::KernelTrait::index 		index;
	typedef typename Problem::KernelTrait::scalar 		scalar;
	typedef typename Problem::KernelTrait::Vector 		Vector;
	typedef typename Problem::KernelTrait::Matrix 		Matrix;

	const Problem& 					__p;
	Matrix 							__mtm;
	Matrix 							__mt;
	scalar							__rho;
	index							__maxiteration;

	/** used variables */
	//%{
	Vector 							__x;
	Vector 							__y;
	Vector 							__z;
	//%}

	/** whether the problem is solved */
	bool 							__is_solved;

#ifdef EIGEN
	Eigen::LDLT<Eigen::MatrixXd>	__linear_solver;
#endif

	
	void doSolve();
	void doCompute();
public:
	/** constructor and deconstructor */
	//%{
	LassoADMMSolver(
		const Problem& p,
		const scalar rho,
		const index maxiteration);
	//%}

	// void compute();
	/** sub-problems */
	//%{
	void xSubProblem(const Vector& zk,const Vector& yk,Vector& xn);
	void zSubProblem(const Vector& xn,const Vector& yk,Vector& zn);
	void ySubProblem(const Vector& xn,const Vector& zn,Vector& yn);
	//%}
	// void solve();

	/** setter and getter */
	//%{
	const Vector& result() const;
	//%}
};

/*		FISTA method for lasso problem
 *
 */
template<class Problem,class Time = NoTimeStatistics>
class FISTALassoSolver
	:
	public GeneralSolver<Time>,
	noncopyable
{
private:
	typedef typename Problem::KernelTrait 			kernel;
	typedef typename kernel::scalar 				scalar;
	typedef typename kernel::index 					index;
	typedef typename kernel::Vector 				Vector;
	typedef typename kernel::Matrix 				Matrix;
	/** private variables */
	//%{
	/** reference to the problem */
	const Problem 			&__p;
	/** Lipschitz constant */
	scalar 					__L;
	/** scaling factor of the constant */
	scalar 					__beta;
	/** whether the operation norm is computed */
	bool 					__optimal_constant;
	/** the result */
	Vector 					__x;
	/** max iteration number */
	index 					__max_iteration;
	/** final iteration number */
	index 					__iter_num;
	/** threshold error */
	scalar 					__thresh;
	/** final error */
	scalar 					__error;
	//%}

	FISTALassoSolver();

public:

	/** constructor and deconstructor */
	//%{
	FISTALassoSolver(
		const Problem& p ,
		const index maxiteration = 10000,
		const scalar thresh = 1e-6,
		const scalar beta = 2.0,
		bool optimalused = false
		);
	//%}

	/** getter */
	//%{
	index finalIterationNumber() const;
	index maxIterationNumber() const;
	scalar threshold() const;
	//%}
};

template<class Problem,class Time>
FISTALassoSolver<Problem,Time>::FISTALassoSolver(
	const Problem& p,
	const index maxiteration,
	const scalar thresh,
	const scalar beta,
	bool optimalused)
	:
	__p(p),
	__beta(beta),
	__optimal_constant(optimalused),
	__max_iteration(maxiteration),
	__iter_num(0),
	__thresh(thresh)
{

}

template<class Problem,class Time>
typename FISTALassoSolver<Problem,Time>::index FISTALassoSolver<Problem,Time>::finalIterationNumber() const
{
	return __iter_num;
}

template<class Problem,class Time>
typename FISTALassoSolver<Problem,Time>::index FISTALassoSolver<Problem,Time>::maxIterationNumber() const
{
	return __max_iteration;
}

template<class Problem,class Time>
typename FISTALassoSolver<Problem,Time>::scalar FISTALassoSolver<Problem,Time>::threshold() const
{
	return __thresh;
}
	
/*		The Lasso problem class
 *
 */
template<class kernel>
class LassoProblem
{
private:
	typedef typename kernel::scalar 		scalar;
	typedef typename kernel::Matrix 		Matrix;
	typedef typename kernel::Vector 		Vector;

	const Matrix& 				__A;
	const Vector& 				__b;
	scalar 						__lambda;
public:
	typedef kernel							KernelTrait;

	LassoProblem( const Matrix& A,
		const Vector& b,
		const scalar lambda = 0.5);

	void proximalSolve();

	/** getter and setter */
	//%{
	const Matrix& matA() const;
	const Vector& obB() const;
	const scalar& lambda() const;
	scalar objective( const Vector& x) const;
	//%}
};


/*********************Implementation of LassoProximalSolver ******************/

template<class Problem,class Time>
void LassoProximalSolver<Problem,Time>::doCompute()
{
	/** compute the beta */
	__s.matA().mtm(__mtm);
	__mt = __s.matA().transpose();
	Eigen::MatrixXd mat(__mtm.rows(),__mtm.cols());
	for ( int i = 0 ; i < __mtm.rows() ; ++ i )
		for ( int j = 0 ; j < __mtm.cols() ; ++ j )
			mat(i,j) = __mtm(i,j);
	__linear_solver.compute(mat);
	scalar e = __s.matA().operationNorm();
	std::cout<<"L is "<<e*e<<std::endl;
	__mu = 1.0/(2*e*e);
}

template<class Problem,class Time>
LassoProximalSolver<Problem,Time>::LassoProximalSolver( 
	const Problem& s , 
	const scalar mu , 
	const scalar beta,
	const index maxiteration)
	:
	__s(s),
	__mu(mu),
	__beta(beta),
	__maxiteration(maxiteration)
{
	doCompute();
}

template<class Problem,class Time>
typename LassoProximalSolver<Problem,Time>::Vector LassoProximalSolver<Problem,Time>::fGradient(const Vector& x)
{
	return 2*(__mtm*x-__mt*__s.obB());
}

template<class Problem,class Time>
typename LassoProximalSolver<Problem,Time>::Vector LassoProximalSolver<Problem,Time>::gProximal( const Vector& x )
{
	Vector result(x.size());
	for ( int i = 0 ; i < x.size() ; ++ i )
	{
		result(i) = computeProximal(AbsFunction(),x(i),__mu*__s.lambda());
	}
	return result;
}

template<class Problem,class Time>
typename LassoProximalSolver<Problem,Time>::scalar LassoProximalSolver<Problem,Time>::f(
	const Vector& x)
{
	return (__s.matA()*x-__s.obB()).squaredNorm();
}

template<class Problem,class Time>
typename LassoProximalSolver<Problem,Time>::scalar LassoProximalSolver<Problem,Time>::fMu(
	const Vector& x,
	const Vector& y)
{
	scalar result = (__s.matA()*y-__s.obB()).squaredNorm();
	result += fGradient(y).dot(x-y);
	result += 1/(2*__mu)*(x-y).squaredNorm();
	return result;
}

template<class Problem,class Time>
bool LassoProximalSolver<Problem,Time>::iteration( const Vector& xk , Vector& xn)
{
	Eigen::VectorXd vec(xk.size());
	Vector temp = fGradient(xk);
	xn = xk-__mu*temp;
	xn = gProximal(xn);
	if (IS_ZERO((xk-xn).squaredNorm()))
	{
		std::cout<<"situation 2 is reached"<<std::endl;
		return true;
	}
	else{
		__mu = __beta*__mu;
		return false;
	}
}

template<class Problem,class Time>
void LassoProximalSolver<Problem,Time>::doSolve()
{
	Vector xp,xn;
	xp = Vector(__s.matA().cols());
	int i = 0;
	while(!iteration(xp,xn)){
		++i;
		if (i>__maxiteration)
			break;
		xp = xn;
	}
	std::cout<<"result is "<<xn<<std::endl;
	std::cout<<i<<" iterations are used!"<<std::endl;
	std::cout<<"final norm is "<<__s.objective(xn)<<std::endl;
}

template<class Problem,class Time>
typename LassoProximalSolver<Problem,Time>::Vector LassoProximalSolver<Problem,Time>::result() const
{
	return __x;
}

/***********************Implementation of ADMM solver*************************/
template<class Problem,class Time>
LassoADMMSolver<Problem,Time>::LassoADMMSolver(
	const Problem& p,
	const scalar rho,
	const index maxiteration)
	:
	__p(p),
	__rho(rho),
	__maxiteration(maxiteration)
{
	this->compute();
}

template<class Problem,class Time>
void LassoADMMSolver<Problem,Time>::doCompute()
{
	__p.matA().mtm(__mtm);
	__mt = __p.matA().transpose();
#ifdef EIGEN
	Eigen::MatrixXd mat(__mtm.rows(),__mtm.cols());
	Eigen::MatrixXd iden = Eigen::MatrixXd::Identity(__mtm.rows(),__mtm.cols());
	for ( int i = 0 ; i < __mtm.rows() ; ++ i )
		iden(i,i) = 1.0/(2*__rho);
	
	for ( int i = 0 ; i < __mtm.rows() ; ++ i )
		for ( int j = 0 ; j < __mtm.cols() ; ++ j )
			mat(i,j) = __mtm(i,j);
	mat = mat+iden;
	__linear_solver.compute(mat);
#endif
}

template<class Problem,class Time>
void LassoADMMSolver<Problem,Time>::xSubProblem(const Vector& zk,const Vector& yk,Vector&xn)
{
	Vector rhd = __mt*__p.obB()+1.0/(2.0*__rho)*(zk-yk);
#ifdef EIGEN
	Eigen::VectorXd vec(rhd.size());
	for ( int i = 0 ; i < rhd.size() ; ++ i )
		vec(i) = rhd(i);
	vec = __linear_solver.solve(vec);
	for ( int i = 0 ; i < rhd.size() ; ++ i )
		xn(i) = vec(i);
#endif
}

template<class Problem,class Time>
void LassoADMMSolver<Problem,Time>::zSubProblem(const Vector& xn,const Vector& yk,Vector& zn)
{
	for ( int i = 0 ; i < zn.size() ; ++ i )
		zn(i) = computeProximal(AbsFunction(),xn(i)+yk(i),__rho*__p.lambda());
}

template<class Problem,class Time>
void LassoADMMSolver<Problem,Time>::ySubProblem(const Vector& xn,const Vector& zn,Vector& yn)
{
	yn = yn+(xn-zn);
}

template<class Problem,class Time>
void LassoADMMSolver<Problem,Time>::doSolve()
{
	index m = __p.matA().rows();
	index n = __p.matA().cols();
	Vector x(n),y(n),z(n),xp(x);
	scalar err;
	int i = 0;
	do{
		xSubProblem(z,y,x);
		zSubProblem(x,y,z);
		ySubProblem(x,z,y);
		err = (x-xp).squaredNorm();
		xp = x;
		if((i++)>=__maxiteration)
			break;
		// std::cout<<err<<std::endl;
	}while(!IS_ZERO(err));
	std::cout<<"result is "<<z<<std::endl;
	std::cout<<i<<" iterations are used "<<std::endl;
	std::cout<<"norm is "<<__p.objective(z)<<std::endl;
}

template<class Problem,class Time>
const typename LassoADMMSolver<Problem,Time>::Vector& LassoADMMSolver<Problem,Time>::result() const
{
	if(!__is_solved)
		std::cerr<<"Warning: the problem is not successfully solved yet!"<<std::endl;
	return __x;
}

/***********************Implementation of LassoProblem************************/

template<class kernel>
LassoProblem<kernel>::LassoProblem(
	const Matrix& A,
	const Vector& b,
	const scalar lambda)
	:
	__A(A),
	__b(b),
	__lambda(lambda)
{  
}

template<class kernel>
void LassoProblem<kernel>::proximalSolve()
{
}

template<class kernel>
const typename kernel::Matrix& LassoProblem<kernel>::matA() const
{
	return __A;
}

template<class kernel>
const typename kernel::Vector& LassoProblem<kernel>::obB() const
{
	return __b;
}

template<class kernel>
const typename kernel::scalar& LassoProblem<kernel>::lambda() const
{
	return __lambda;
}

template<class kernel>
typename kernel::scalar LassoProblem<kernel>::objective( const Vector& x ) const
{
	return (__A*x-__b).squaredNorm()+__lambda*x.absNorm();
}
}

#endif