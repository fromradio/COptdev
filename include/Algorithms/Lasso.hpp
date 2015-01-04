// This file is part of COPT, a light-weight C++ based optimization open source library
//
// Copyright (C) 2015 Ruimin Wang <ruimin.wang13@gmail.com>
// Copyright (C) 2015 MathU
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#ifndef LASSO_HPP__
#define LASSO_HPP__

namespace COPT
{
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
	const Problem& 								__p;
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
	//%}

	// void init();

	LassoProximalSolver();

	void solvingBegin();
	void doCompute();
	scalar doOneIteration();
	

public:


	typedef proximal_solver 				ObjectCategory;
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
	
	/** setter and getter */
	//%{
	/** the final result of proximal solver */
	const Vector& result() const;
	//%}

	scalar objective() const;

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
	// index							__maxiteration;

	/** used variables */
	//%{
	Vector 							__x;
	Vector 							__y;
	Vector 							__z;
	//%}

	/** whether the problem is solved */
	bool 							__is_solved;

	/** linear solver */
	LU<Matrix> 						__linear_solver;

	/** implementation of virtual functions */
	//%{
	void solvingBegin();
	void doCompute();
	scalar doOneIteration();
	scalar objective() const;
	//%}

public:

	typedef admm_solver 				ObjectCategory;
	/** constructor and deconstructor */
	//%{
	LassoADMMSolver(
		const Problem& p,
		const scalar rho = 0.5,
		const index maxiteration = 1000,
		const scalar thresh = 1e-8);
	//%}

	/** sub-problems */
	//%{
	void xSubProblem(const Vector& zk,const Vector& yk,Vector& xn);
	void zSubProblem(const Vector& xn,const Vector& yk,Vector& zn);
	void ySubProblem(const Vector& xn,const Vector& zn,Vector& yn);
	//%}

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
	/** parameter t */
	scalar 					__t;

	FISTALassoSolver();


	void solvingBegin();
	void doCompute();
	scalar doOneIteration();
	scalar objective() const;

public:

	typedef fista_solver 					ObjectCategory;
	/** constructor and deconstructor */
	//%{
	FISTALassoSolver(
		const Problem& p ,
		const scalar L = 1.0,
		const index maxiteration = 100000,
		const scalar thresh = 1e-6,
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
	const Vector& result() const;
	//%}
};


template <class kernel>
class VectorProblem
{

	typedef typename kernel::scalar 		scalar;
	typedef typename kernel::podscalar 		podscalar;
	typedef typename kernel::index 			index;
	typedef typename kernel::Matrix 		Matrix;
	typedef typename kernel::Vector 		Vector;

	/** the dimension of the problem */
	index 					__dim;

public:

	typedef kernel 							KernelTrait;

	/** constructor and deconstructor */
	//%{
	/** default constructor */
	VectorProblem( const index dim = 0 );
	virtual ~VectorProblem(){}
	//%}

	/** compute the objective function */
	virtual podscalar objective( const Vector& x ) const = 0;

	/** check whether the input is valid */
	virtual bool isValidInput( const Vector& x ) const = 0;

	/** check whether the problem is a valid problem */
	virtual bool isValid( ) const = 0;

	/** validation of the problem */
	virtual void validation() const;

	/** getter and setter */
	//%{
	void setDimension( const index dim );
	index dimension() const;
	//%}
};

template<class kernel>
VectorProblem<kernel>::VectorProblem( const index dim )
	:
	__dim(dim)
{
}

template<class kernel>
void VectorProblem<kernel>::validation() const
{
	if(!this->isValid())
	{
		throw COException("Vector problem error: Validation of this problem does not pass!");
	}
}

template<class kernel>
void VectorProblem<kernel>::setDimension( const index dim )
{
	__dim = dim;
}

template<class kernel>
typename VectorProblem<kernel>::index VectorProblem<kernel>::dimension() const
{
	return __dim;
}
	
/*		The Lasso problem class
 *
 */
template<class kernel>
class LassoProblem
	:
	public VectorProblem<kernel>
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
	typedef lasso_problem					ObjectCategory;

	LassoProblem(
		const Matrix& A,
		const Vector& b,
		const scalar lambda = 0.5);

	void proximalSolve();

	/** check whether the input is valid */
	bool isValidInput( const Vector& x ) const;

	/** check whether the problem is a valid problem */
	bool isValid( ) const;

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
	__p.matA().mtm(__mtm);
	__mt = __p.matA().transpose();
	scalar e = __p.matA().operationNorm();
	__mu = 1.0/(2*e*e);
}

template<class Problem,class Time>
LassoProximalSolver<Problem,Time>::LassoProximalSolver( 
	const Problem& s , 
	const scalar mu , 
	const scalar beta,
	const index maxiteration)
	:
	__p(s),
	__mu(mu),
	__beta(beta),
	__maxiteration(maxiteration)
{
	this->doCompute();
}

template<class Problem,class Time>
typename LassoProximalSolver<Problem,Time>::Vector LassoProximalSolver<Problem,Time>::fGradient(const Vector& x)
{
	return 2*(__mtm*x-__mt*__p.obB());
}

template<class Problem,class Time>
typename LassoProximalSolver<Problem,Time>::Vector LassoProximalSolver<Problem,Time>::gProximal( const Vector& x )
{
	Vector result(x.size());
	for ( int i = 0 ; i < x.size() ; ++ i )
	{
		result(i) = computeProximal(AbsFunction(),x(i),__mu*__p.lambda());
	}
	return result;
}

template<class Problem,class Time>
typename LassoProximalSolver<Problem,Time>::scalar LassoProximalSolver<Problem,Time>::f(
	const Vector& x)
{
	return (__p.matA()*x-__p.obB()).squaredNorm();
}

template<class Problem,class Time>
typename LassoProximalSolver<Problem,Time>::scalar LassoProximalSolver<Problem,Time>::fMu(
	const Vector& x,
	const Vector& y)
{
	scalar result = (__p.matA()*y-__p.obB()).squaredNorm();
	result += fGradient(y).dot(x-y);
	result += 1/(2*__mu)*(x-y).squaredNorm();
	return result;
}

template<class Problem,class Time>
void LassoProximalSolver<Problem,Time>::solvingBegin()
{
	__x.resize(__p.matA().cols());
}
template<class Problem,class Time>
typename LassoProximalSolver<Problem,Time>::scalar LassoProximalSolver<Problem,Time>::objective() const
{
	return __p.objective(__x);
}

template<class Problem,class Time>
typename LassoProximalSolver<Problem,Time>::scalar LassoProximalSolver<Problem,Time>::doOneIteration()
{
	// Eigen::VectorXd vec(xk.size());
	Vector xp = __x;
	Vector temp = fGradient(__x);
	__x = __x-__mu*temp;
	__x = gProximal(__x);
	return std::sqrt((xp-__x).squaredNorm());
}

template<class Problem,class Time>
const typename LassoProximalSolver<Problem,Time>::Vector& LassoProximalSolver<Problem,Time>::result() const
{
	return __x;
}

/***********************Implementation of ADMM solver*************************/
template<class Problem,class Time>
LassoADMMSolver<Problem,Time>::LassoADMMSolver(
	const Problem& p,
	const scalar rho,
	const index maxiteration,
	const scalar thresh)
	:
	GeneralSolver<typename Problem::KernelTrait,Time>(maxiteration,thresh),
	__p(p),
	__rho(rho)
{
	this->compute();
}

template<class Problem,class Time>
void LassoADMMSolver<Problem,Time>::doCompute()
{
	__p.matA().mtm(__mtm);
	__mt = __p.matA().transpose();
	Matrix mat(__mtm);
	for ( int i = 0 ; i < mat.rows() ; ++ i )
	{
		mat(i,i) += 1.0/(2*__rho);	
	}
	__linear_solver.compute(mat);
}

template<class Problem,class Time>
void LassoADMMSolver<Problem,Time>::xSubProblem(const Vector& zk,const Vector& yk,Vector&xn)
{
	Vector rhd = __mt*__p.obB()+1.0/(2.0*__rho)*(zk-yk);
	xn = __linear_solver.solve(rhd);
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
void LassoADMMSolver<Problem,Time>::solvingBegin()
{
	__x.resize(__p.matA().cols());
	__y.resize(__p.matA().cols());
	__z.resize(__p.matA().cols());
}

template<class Problem,class Time>
typename LassoADMMSolver<Problem,Time>::scalar LassoADMMSolver<Problem,Time>::doOneIteration()
{
	Vector xp=__x;
	xSubProblem(__z,__y,__x);
	zSubProblem(__x,__y,__z);
	ySubProblem(__x,__z,__y);
	return std::abs(__p.objective(__x)-__p.objective(xp));
}

template<class Problem,class Time>
typename LassoADMMSolver<Problem,Time>::scalar LassoADMMSolver<Problem,Time>::objective() const
{
	return __p.objective(__x);
}

template<class Problem,class Time>
const typename LassoADMMSolver<Problem,Time>::Vector& LassoADMMSolver<Problem,Time>::result() const
{
	if(!__is_solved)
		std::cerr<<"Warning: the problem is not successfully solved yet!"<<std::endl;
	return __x;
}

/**********************Implementation of FISTALASSOSolver*********************/

template<class Problem,class Time>
void FISTALassoSolver<Problem,Time>::doCompute()
{
	__x.resize(__p.matA().cols());
	__y.resize(__p.matA().cols());
	__atb = __p.matA().transMulti(__p.obB());
	scalar e = __p.matA().operationNorm();
	__L = 2*e*e;
	__optimal_constant = true;
}

template<class Problem,class Time>
void FISTALassoSolver<Problem,Time>::solvingBegin()
{
	__t = 1.0;
}

template<class Problem,class Time>
typename FISTALassoSolver<Problem,Time>::scalar FISTALassoSolver<Problem,Time>::doOneIteration()
{
	Vector ply,xp = __x;
	pArgmin(__y,__L,ply);
	if(!__optimal_constant)
		__L = findConstant(ply,__y,__L,__beta);
	pArgmin(__y,__L,__x);
	scalar tn =  (1.0+std::sqrt(1+4*__t*__t))/2.0;
	__y = __x +((__t-1)/tn)*(__x-xp);
	__t = tn;
	return std::abs(__p.objective(__x)-__p.objective(xp));
}

template<class Problem,class Time>
typename FISTALassoSolver<Problem,Time>::scalar FISTALassoSolver<Problem,Time>::objective() const
{
	return __p.objective(__x);
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
	GeneralSolver<kernel,Time>(maxiteration,thresh),
	__p(p),
	__L(L),
	__beta(beta),
	__optimal_constant(optimalused)
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

template<class Problem,class Time>
const typename FISTALassoSolver<Problem,Time>::Vector& FISTALassoSolver<Problem,Time>::result() const
{
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
bool LassoProblem<kernel>::isValid( ) const
{
	if (__A.rows()!=__b.size())
		return false;
	else 
		return true;
}

template<class kernel>
bool LassoProblem<kernel>::isValidInput( const Vector& x ) const
{
	if (__A.cols() != x.size() )
		return false;
	else 
		return true;
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