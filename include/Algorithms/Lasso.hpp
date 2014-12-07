//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef LASSO_HPP__
#define LASSO_HPP__

namespace COPT
{

template<class Problem>
class LassoProximalSolver
{
private:
	typedef typename Problem::KernelTrait::index		index;
	typedef typename Problem::KernelTrait::scalar 		scalar;
	typedef typename Problem::KernelTrait::Vector 		Vector;
	typedef typename Problem::KernelTrait::Matrix 		Matrix;

	/** the corresponding lasso problem */
	const Problem& 				__s;
	Matrix 						__mtm;
	Matrix 						__mt;
	scalar 						__mu;
	scalar 						__beta;
	index 						__maxiteration;

	Eigen::LDLT<Eigen::MatrixXd>				__linear_solver;

	void init();
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

	/** solve the problem */
	void solve();

};

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

template<class Problem>
void LassoProximalSolver<Problem>::init()
{
	__s.matA().mtm(__mtm);
	__mt = __s.matA().transpose();
	Eigen::MatrixXd mat(__mtm.rows(),__mtm.cols());
	for ( int i = 0 ; i < __mtm.rows() ; ++ i )
		for ( int j = 0 ; j < __mtm.cols() ; ++ j )
			mat(i,j) = __mtm(i,j);
	__linear_solver.compute(mat);
}

template<class Problem>
LassoProximalSolver<Problem>::LassoProximalSolver( 
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
	init();
}

template<class Problem>
typename LassoProximalSolver<Problem>::Vector LassoProximalSolver<Problem>::fGradient(const Vector& x)
{
	return __mtm*x-__mt*__s.obB();
}

template<class Problem>
typename LassoProximalSolver<Problem>::Vector LassoProximalSolver<Problem>::gProximal( const Vector& x )
{
	Vector result(x.size());
	for ( int i = 0 ; i < x.size() ; ++ i )
	{
		result(i) = computeProximal(AbsFunction(),x(i),__mu*__s.lambda());
	}
	return result;
}

template<class Problem>
typename LassoProximalSolver<Problem>::scalar LassoProximalSolver<Problem>::f(
	const Vector& x)
{
	return 0.5*(__s.matA()*x-__s.obB()).squaredNorm();
}

template<class Problem>
typename LassoProximalSolver<Problem>::scalar LassoProximalSolver<Problem>::fMu(
	const Vector& x,
	const Vector& y)
{
	scalar result = 0.5*(__s.matA()*y-__s.obB()).squaredNorm();
	result += fGradient(y).dot(x-y);
	result += 1/(2*__mu)*(x-y).squaredNorm();
	return result;
}

template<class Problem>
bool LassoProximalSolver<Problem>::iteration( const Vector& xk , Vector& xn)
{
	// std::cout<<"gradient is "<<fGradient(xk)<<std::endl;
	Eigen::VectorXd vec(xk.size());
	Vector temp = fGradient(xk);
	for ( int i = 0 ; i < temp.size() ;++ i )
		vec(i) = temp(i);
	vec = __linear_solver.solve(vec);
	for ( int i = 0 ; i < temp.size() ; ++ i )
		temp(i) = vec(i);
	xn = xk-__mu*temp;
	// xn = gProximal(xk-__mu*fGradient(xk));
	// std::cout<<xn<<std::endl;
	scalar f1 = f(xn) , f2 = fMu(xn,xk);
	if ( f1 <= f2 )
	{
		std::cout<<"situation 1 is reached!"<<std::endl;
		return true;
	}
	else if (IS_ZERO(f1-f2))
	{
		std::cout<<"situation 2 is reached"<<std::endl;
		return true;
	}
	else{
		__mu = __beta*__mu;
		return false;
	}
}

template<class Problem>
void LassoProximalSolver<Problem>::solve()
{
	Vector xp,xn;
	xp = Vector(__s.matA().cols());
	int i = 0;
	while(!iteration(xp,xn)){
		++i;
		if (i>__maxiteration)
			break;
		xp = xn;
		std::cout<<i<<std::endl;
	}
	std::cout<<"result is "<<xn<<std::endl;
	std::cout<<i<<" iterations are used!"<<std::endl;
	std::cout<<"final norm is "<<__s.objective(xn)<<std::endl;
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
	return 0.5*(__A*x-__b).squaredNorm()+__lambda*x.absNorm();
}
}

#endif