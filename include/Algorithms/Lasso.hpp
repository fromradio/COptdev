//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef LASSO_HPP__
#define LASSO_HPP__

namespace COPT
{

/*		A general design for solver. The solver derives from a Time
 *		Stastistics class to help it compute time cost of the solver.
 *		The solver contains general information like max iteration 
 *		number, final iteration number, threshold, estimated error
 */	
template<class kernel , class Time=NoTimeStatistics>
class GeneralSolver
	:
	public COPTObject,
	public Time,
	public noncopyable
{
	typedef typename kernel::index 				index;
	typedef typename kernel::scalar 			scalar;

	/** max iteration number */
	index 		__max_iteration;
	/** final iteration number */
	index 		__iter_num;
	/** error threshold */
	scalar 		__thresh;
	/** final estimated error */
	scalar 		__estimated_error;

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

	/** getter and setter */
	//%{
	void 	setMaxIteration(const index maxiteration);
	index 	maxIterationNumber() const;
	void 	setIterationNumber(const index iter);
	index 	iterationNumber() const;
	void 	setThreshold( const scalar thresh);
	scalar 	threshold() const;
	void 	setEstimatedError(const scalar erro);
	scalar 	estimatedError() const;
	//%}
};

template<class kernel,class Time>
void GeneralSolver<kernel,Time>::compute()
{
	this->computationBegin();
	this->doCompute();
	this->computationEnd();
}

template<class kernel,class Time>
void GeneralSolver<kernel,Time>::solve()
{
	this->solvingBegin();
	this->doSolve();
	this->solvingEnd();
	this->printTimeInfo();
}

template<class kernel,class Time>
void GeneralSolver<kernel,Time>::setMaxIteration( const index maxiteration)
{
	__max_iteration = maxiteration;
}

template<class kernel,class Time>
typename GeneralSolver<kernel,Time>::index GeneralSolver<kernel,Time>::maxIterationNumber() const
{
	return __max_iteration;
}

template<class kernel,class Time>
void GeneralSolver<kernel,Time>::setIterationNumber( const index iter )
{
	__iter_num = iter;
}

template<class kernel,class Time>
typename GeneralSolver<kernel,Time>::index GeneralSolver<kernel,Time>::iterationNumber() const
{
	return __iter_num;
}

template<class kernel,class Time>
void GeneralSolver<kernel,Time>::setThreshold( const scalar thresh )
{
	__thresh = thresh;
}

template<class kernel,class Time>
typename GeneralSolver<kernel,Time>::scalar GeneralSolver<kernel,Time>::threshold() const
{
	return __thresh;
}

template<class kernel,class Time>
void GeneralSolver<kernel,Time>::setEstimatedError(const scalar error )
{
	__estimated_error = error;
}

template<class kernel,class Time>
typename GeneralSolver<kernel,Time>::scalar GeneralSolver<kernel,Time>::estimatedError() const
{
	return __estimated_error;
}


/*		Proximal solver for lasso problem
 *
 *
 */
template<class Problem,class Time = NoTimeStatistics>
class LassoProximalSolver
	:
	public GeneralSolver<typename Problem::KernelTrait,Time>
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
	public GeneralSolver<typename Problem::KernelTrait,Time>
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
	public GeneralSolver<typename Problem::KernelTrait,Time>
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
	/** the y vector */
	Vector 					__y;
	/** A^Tb vector */
	Vector 					__atb;
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

	void doSolve();
	void doCompute();

public:

	/** constructor and deconstructor */
	//%{
	FISTALassoSolver(
		const Problem& p ,
		const scalar L = 1.0,
		const index maxiteration = 100000,
		const scalar thresh = 1e-12,
		const scalar beta = 2.0,
		bool optimalused = false
		);
	//%}

	inline scalar fFunction( const Vector& x );
	inline Vector fGradient( const Vector& x );
	inline scalar qFunction( const Vector& x , const Vector& y , const scalar L);

	/** compute the p_L(y) problem of solver */
	inline void pArgmin( const Vector& y , const scalar L , Vector& x);
	/** find the next L */
	inline scalar findConstant(const Vector& ply , const Vector& y , const scalar L , const scalar beta );

	/** getter */
	//%{
	//%}
};

template<class Problem,class Time>
void FISTALassoSolver<Problem,Time>::doCompute()
{
	std::cout<<"computation"<<std::endl;
	__x.resize(__p.matA().cols());
	__y.resize(__p.matA().cols());
	__atb = __p.matA().transMulti(__p.obB());
	scalar e = __p.matA().operationNorm();
	__L = 2*e*e;
	__optimal_constant = true;
}

template<class Problem,class Time>
void FISTALassoSolver<Problem,Time>::doSolve()
{
	scalar t = 1.0;
	
	Vector ply;
	do
	{
		Vector xp=__x;
		pArgmin(__y,__L,ply);
		if(!__optimal_constant)
			__L = findConstant(ply,__y,__L,__beta);
		pArgmin(__y,__L,__x);
		scalar tn = (1.0+std::sqrt(1+4*t*t))/2.0;
		__y = __x+((t-1)/tn)*(__x-xp);
		t = tn;
		__error = (__x-xp).squaredNorm();
		if (__iter_num ++ >= __max_iteration )
			break;
	} while(__error>__thresh);
	std::cout<<"iteration number is "<<__iter_num<<std::endl;
	std::cout<<"result is "<<__x<<std::endl;
	std::cout<<"norm is "<<__p.objective(__x)<<std::endl;
	std::cout<<"final L is "<<__L<<std::endl;
}

template<class Problem,class Time>
FISTALassoSolver<Problem,Time>::FISTALassoSolver(
	const Problem& p,
	const scalar L,
	const index maxiteration,
	const scalar thresh,
	const scalar beta,
	bool optimalused)
	:
	__p(p),
	__L(L),
	__beta(beta),
	__optimal_constant(optimalused),
	__max_iteration(maxiteration),
	__iter_num(0),
	__thresh(thresh)
{
	this->compute();
}

template<class Problem,class Time>
typename FISTALassoSolver<Problem,Time>::scalar FISTALassoSolver<Problem,Time>::fFunction( const Vector& x )
{
	return (__p.matA()*x-__p.obB()).squaredNorm();
}

template<class Problem,class Time>
typename FISTALassoSolver<Problem,Time>::Vector FISTALassoSolver<Problem,Time>::fGradient( const Vector & x )
{
	return 2*(__p.matA().transMulti(__p.matA()*x)-__atb);
}

template<class Problem,class Time>
typename FISTALassoSolver<Problem,Time>::scalar FISTALassoSolver<Problem,Time>::qFunction( const Vector& x , const Vector& y , const scalar L )
{
	scalar q = fFunction(y);
	q += (x-y).dot(fGradient(y))+0.5*L*(x-y).squaredNorm();
	return q;
}

template<class Problem,class Time>
void FISTALassoSolver<Problem,Time>::pArgmin( const Vector& y , const scalar L , Vector& x )
{
	x.resize(y.size());
	Vector v = y-1.0/L*fGradient(y);
	scalar factor = __p.lambda()/L;
	for ( int i = 0 ; i < x.size() ; ++ i )
	{
		x(i) = computeProximal(AbsFunction(),v(i),factor);
	}
}

template<class Problem,class Time>
typename FISTALassoSolver<Problem,Time>::scalar FISTALassoSolver<Problem,Time>::findConstant( const Vector& ply , const Vector& y , const scalar L , const scalar beta )
{
	if ( beta < 1.0 )
		throw COException("please make sure the ratio factor beta of FISTA solver is bigger than 1.0!");
	scalar nL = L;
	while ( fFunction(ply) > qFunction(ply,y,nL) )
	{
		nL = beta*nL;
	}
	return nL;
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
	this->compute();
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