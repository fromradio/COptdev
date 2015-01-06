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


#ifndef SIMPLEX_HPP__
#define SIMPLEX_HPP__


namespace COPT
{


/*		The Linear Programming(LP) problem is assuemd to be with the
 *		standard form. The constraint contains two parts: Ax=b and x>=0.
 *		Any non-standard problem should be transformed at first.
 *		The very basic algorithm is implemented. We use a two phase 
 *		algorithm. The first phase is used to find a feasible x.
 *		A linear programming problem takes weight vector, linear constraint
 *		as input. Output information contains final status of solver, success,
 *		not feasible or unbounded.
 */
template<class kernel>
class SimplexSolver
{
private:
	typedef typename kernel::scalar				scalar;
	typedef typename kernel::index 				index;
	typedef typename kernel::Vector				Vector;
	typedef typename kernel::Matrix				Matrix;

	/** type of the result of a step */
	enum StepType{
		Normal,
		Optimal,
		Unbounded
	};

	/** type of the algorithm */
	enum AlgoType{
		Starting,
		Success,
		Unbound,
		NotFeasible
	};

	/** the variable of the problem*/
	//%{
	/** the dimension of the problem */
	index 				__n;
	/** the dimension of right hand vector b */
	index 				__m;
	/** the matrix A in equal constraint */
	Matrix 				__A;
	/** the right hand vector b */
	Vector 				__b;
	/** the weight vector in target function */
	Vector 				__c;
	/** the variable */
	Vector				__x;
	/** the final evaluation of the function */
	scalar				__val;
	/** the status of the solver */
	AlgoType 			__type;
	//%} end of variables

public:
	/** Constructor and Deconstructor */
	//%{
	/** default constructor */
	SimplexSolver();
	/** standard constructor */
	SimplexSolver(
		const Matrix& A,
		const Vector& b,
		const Vector& c);
	//%}

	/** judge the consistency of the problem */
	bool consistencyJudgement() const;

	/** solve the linear programming problem */
	AlgoType solve();

	/** getter and setter of SimplexSolver */
	//%{
	/** the current status solver */
	AlgoType solverStatus() const;
	/** a string containing the current status of the solver */
	std::string stringSolverStatus() const;
	/** matrix A */
	const Matrix& matA() const;
	/** right hand vector b */
	const Vector& rhb() const;
	/** weight vector c */
	const Vector& weight() const;
	/** result x */
	const Vector& result() const;
	/** evaluation value */
	const scalar& val() const;
	//%} end of getter and setter

public:

	/** static inline functions of simplex solver */
	//%{
	/* 		Find the solution according to the matrix and blocking vector
	 *		/param A:			input matrix
	 *		/param b:			right hand vector
	 *		/param indb:		the blocking vector
	 *		/param x:			the result on output
	 */
	static inline void solveBlockingSystem(
						const Matrix& A,
						const Vector& b,
						const std::vector<index>& indb,
						Vector& x);

	/*		Find the complete solution with a linear system and blocking vector
	 *		/param A:			input matrix
	 *		/param b:			right hand vector
	 *		/param indb:		the blocking vector
	 *		/param x:			the result on output
	 */
	static inline void solveCompleteBlockingSystem(
						const Matrix& A,
						const Vector& b,
						const std::vector<index>& indb,
						Vector& x);

	/*		One step of a simplex algorithm for solving a Linear Programming
 	 *		problem.
 	 *		/param A: 			input matrix
 	 *		/param b: 			right hand vector
 	 *		/param c: 			weight matrix
 	 *		/param indb:		index of B
 	 *		/param indn:		index of N
 	 */
	static inline StepType oneSimplexStep(
 						const Matrix& A,
 						const Vector& b,
 						const Vector& c,
 						std::vector<index>& indb,
 						std::vector<index>& indn);

 	static inline bool findPivoting(
 						const Vector& xb,
 						const Vector& d,
 						const int q,
 						std::vector<index>& indb,
 						std::vector<index>& indn);

 	/* 		The first phase of simplex method. The target is to find the feasible
	 *		point of the problem. 
	 *				min e^Tz
	 *				s.t. Ax+Ez=b
	 *					 (x,z)>=0
	 *		/param m:
	 *		/param n:
	 *		/param A:
	 *		/param b:
	 *		/param indb:
	 *		/param indn:
	 *		/param indz:
	 */
	static inline bool findFeasiblePoint(
						const index m,
						const index n,
						const Matrix& A,
						const Vector& b,
						std::vector<index>& indb,
						std::vector<index>& indn,
						std::vector<index>& indz);


 	/* 		Simplex method algorithm when initial point is known.
 	 *		/param A:			equal constraint matrix
 	 *		/param b:			right hand vector of equal constraint
 	 *		/param c:			weight vector
 	 *		/param x:			initial feasible point on input and result on output
 	 *		/param val:			the evaluation of the problem
 	 *		/param indb:		initial B on input and final B on output
 	 *		/param indn:		initial N on input and final N on output
 	 */
 	static inline AlgoType simplexMethodWithInitalization(
 						const int m,
 						const int n,
 						const Matrix& A,
 						const Vector& b,
 						const Vector& c,
 						Vector& x,
 						scalar& val,
 						std::vector<index>& indb,
 						std::vector<index>& indn);

 	/*		Phase two of simplex method:
 	 *		The input contains the current indices of used B. indz indicates
 	 *		indices of artificial variable z.
 	 *		/param m:			the dimension of b
 	 *		/param n:			the dimension of x
 	 *		/param A:			input matrix of equal constraint
 	 *		/param b:			right hand vector of equal constraint
 	 *		/param c:			weight vector
 	 *		/param x:			the result on output
 	 *		/param val:			evaluation of the problem
 	 *		/param indb:		currend B
 	 *		/param indn:		current N
 	 *		/param indz:		used Z
 	 */
 	static inline AlgoType simplexMethodWithIndices(
 						const int m,
 						const int n,
 						const Matrix& A,
 						const Vector& b,
 						const Vector& c,
 						Vector& x,
 						scalar& val,
 						std::vector<index>& indb,
 						std::vector<index>& indn,
 						std::vector<index>& indz);
 	//%} end of static methods
}; // End of class SimplexSolver




/** implementation part of SimplexSolver*/

template<class kernel>
SimplexSolver<kernel>::SimplexSolver()
	:
	__m(0),
	__n(0),
	__type(Starting)
{
}

template<class kernel>
SimplexSolver<kernel>::SimplexSolver(
	const Matrix& A,
	const Vector& b,
	const Vector& c)
	:
	__n(A.cols()),
	__m(A.rows()),
	__A(A),
	__b(b),
	__c(c),
	__type(Starting)
{
	if(!consistencyJudgement())
		throw COException("The input of linear programming is not consistent");
}

/** getter and setter */
template<class kernel>
typename SimplexSolver<kernel>::AlgoType SimplexSolver<kernel>::solverStatus() const
{
	return __type;
}

template<class kernel>
std::string SimplexSolver<kernel>::stringSolverStatus() const
{
	switch(__type){
		case Starting:
		{
			return std::string("algorithm has not started yet.");
		}
		break;
		case Success:
		{
			return std::string("algorithm success.");
		}
		break;
		case NotFeasible:
		{
			return std::string("there is no feasible point.");
		}
		break;
		case Unbound:
		{
			return std::string("there is no bound of the problem.");
		}
		break;
		default:
		{
			return std::string("Unknown type");
		}
		break;
	}
}

template<class kernel>
const typename kernel::Matrix& SimplexSolver<kernel>::matA() const
{
	return __A;
}

template<class kernel>
const typename kernel::Vector& SimplexSolver<kernel>::rhb() const
{
	return __b;
}

template<class kernel>
const typename kernel::Vector& SimplexSolver<kernel>::weight() const
{
	return __c;
}

template<class kernel>
const typename kernel::Vector& SimplexSolver<kernel>::result() const
{
	return __x;
}

template<class kernel>
const typename kernel::scalar& SimplexSolver<kernel>::val() const
{
	return __val;
}

template<class kernel>
bool SimplexSolver<kernel>::consistencyJudgement() const
{
	if(__m!=__A.rows()) return false;
	if(__m!=__b.size()) return false;
	if(__n!=__A.cols()) return false;
	if(__n!=__c.size()) return false;
	return true;
}

template<class kernel>
typename SimplexSolver<kernel>::AlgoType SimplexSolver<kernel>::solve()
{
	std::vector<index> indb,indn,indz;
	if( findFeasiblePoint(__m,__n,__A,__b,indb,indn,indz) )
	{
		// the problem is feasible
		__type = simplexMethodWithIndices(__m,__n,__A,__b,__c,__x,__val,indb,indn,indz);
		return __type;
	}
	else{
		__type = NotFeasible;
		return NotFeasible;
	}
}

template<class kernel>
void SimplexSolver<kernel>::solveBlockingSystem(
	const Matrix& A,
	const Vector& b,
	const std::vector<index>& indb,
	Vector& x)
{
	Matrix B;
	B.columnBlockFromMatrix(A,indb);
	x = B.solve(b);
}

template<class kernel>
void SimplexSolver<kernel>::solveCompleteBlockingSystem(
	const Matrix& A,
	const Vector& b,
	const std::vector<index>& indb,
	Vector& x)
{
	x.resize(A.cols());
	Vector xb;
	solveBlockingSystem(A,b,indb,xb);
	for ( int i = 0 ; i < indb.size() ; ++ i )
		x(indb[i]) = xb(i);
}

template<class kernel>
bool SimplexSolver<kernel>::findPivoting(const Vector& xb,const Vector& d,const int q,std::vector<index>& indb,std::vector<index>& indn)
{
	int p = -1;
	scalar min = INFTY;
	for ( int i = 0 ; i < d.size() ; ++ i ){
		if( d[i] > 0 ){
			if(xb[i]/d[i]<min){
				p = i;
				min = xb[i]/d[i];
			}
		}
	}
	if (p == -1 ){
		return false;
	}
	else{
		std::swap(indb[p],indn[q]);
		return true;
	}
}

template<class kernel>
typename SimplexSolver<kernel>::StepType SimplexSolver<kernel>::oneSimplexStep(
	const Matrix& A,
	const Vector& b,
	const Vector& c,
	std::vector<index>& indb,
	std::vector<index>& indn
	)
{
	Matrix B,N;
	B.columnBlockFromMatrix(A,indb);
	Vector xb = B.solve(b);
	Vector lambda = (B.transpose()).solve(c.block(indb));
	Vector cn = c.block(indn);
	N.columnBlockFromMatrix(A,indn);
	Vector sn = cn-N.transpose()*lambda;
	int q = -1;
	for ( int i = 0 ; i < sn.size() ; ++ i ){
		if(sn[i]<0){
			q = i;
			break;
		}
	} // find the corresponding q
	if ( q == -1 )
		return Optimal; // optimal solution is found
	Vector Aq = A.col(indn[q]);
	Vector d = B.solve(Aq);
	bool unbound = true;
	for ( int i = 0 ; i < d.size() ; ++ i ){
		if ( d[i] > 0 )
			unbound = false;
	}
	if (unbound)
		return Unbounded;
	if(!findPivoting(xb,d,q,indb,indn))
		throw COException("Unknown error in simplex solver!");
	return Normal;
}


template<class kernel>
bool SimplexSolver<kernel>::findFeasiblePoint(
	const index m,
	const index n,
	const Matrix& A,
	const Vector& b,
	std::vector<index>& indb,
	std::vector<index>& indn,
	std::vector<index>& indz)
{
	Matrix E = Matrix::identity(m,m),AE;
	for ( int i = 0 ; i < m ; ++ i ){
		if(b(i)<0)
			E(i,i) = static_cast<scalar>(-1.0);
		else
			E(i,i) = static_cast<scalar>(1.0);
	}
	Matrix::stCombineAlongColumn(A,E,AE);
	Vector e(m+n);
	for ( int i = 0 ; i < m ; ++ i )
		e(i+n)=1.0;
	Vector xz(m+n);
	for ( int i = 0 ; i < m ; ++ i )
		xz(i+n) = std::abs(b(i));
	std::vector<index> iindb(m),iindn(n); // the indices for B and N
	for ( index i = 0 ; i < m ; ++ i )
		iindb[i] = n+i;
	for ( index i = 0 ; i < n ; ++ i )
		iindn[i] = i;
	Vector x;
	scalar val;
	AlgoType t = simplexMethodWithInitalization(m,n,AE,b,e,x,val,iindb,iindn);

	// compute indb and indn. index less than n is kept
	indb.clear();indb.reserve(iindb.size());
	indn.clear();indn.reserve(iindn.size());
	indz.clear();indz.reserve(iindb.size());
	for ( index i = 0 ; i < iindb.size() ; ++ i )
	{
		if(iindb[i]<n)
			indb.push_back(iindb[i]);
		else
			indz.push_back(iindb[i]-n);
	}
	for (index i = 0 ; i < iindn.size() ; ++ i )
	{
		if(iindn[i]<n)
			indn.push_back(iindn[i]);
	}

	if ( t == Unbound || t == NotFeasible )
		return false;
	else if ( IS_ZERO(val) )
		return true;
	else
		return false;
}

template<class kernel>
typename SimplexSolver<kernel>::AlgoType SimplexSolver<kernel>::simplexMethodWithInitalization(
	const int m,
	const int n,
	const Matrix& A,
	const Vector& b,
	const Vector& c,
	Vector& x,
	scalar& val,
	std::vector<index>&indb,
	std::vector<index>&indn)
{
	StepType type;
	do{
		type = oneSimplexStep(A,b,c,indb,indn);
	}while(type == Normal);
	if ( type == Unbounded)
		return Unbound;
	else{
		solveCompleteBlockingSystem(A,b,indb,x);
		val = x.dot(c);
		return Success;
	}
}

template<class kernel>
typename SimplexSolver<kernel>::AlgoType SimplexSolver<kernel>::simplexMethodWithIndices(
	const int m,
	const int n,
	const Matrix& A,
	const Vector& b,
	const Vector& c,
	Vector& x,
	scalar& val,
	std::vector<index>&indb,
	std::vector<index>&indn,
	std::vector<index>&indz)
{
	// the difference is that some artificial variables still remains
	index zsize = indz.size();
	// new c
	Vector nc(n+2*zsize);
	for ( int i = 0 ; i < n ; ++ i )
		nc(i)=c(i);
	Matrix nA(m+zsize,n+2*zsize);
	for ( int i = 0 ; i < m ; ++ i )
		for ( int j = 0 ; j < n ; ++ j )
			nA(i,j) = A(i,j);
	for ( int i = 0 ; i < zsize ; ++ i )
	{
		nA(indz[i],n+i) = 1.0;
		nA(m+i,n+i) = 1.0;
		nA(m+i,n+zsize+i) = 1.0;
	}
	// new b
	Vector nb(m+zsize);
	for ( int i = 0 ; i < m ; ++ i )
		nb(i) = b(i);
	Vector nx;
	std::vector<index> iindb(indb);
	for ( int i = 0 ; i < zsize ; ++ i )
		iindb.push_back(n+i);
	AlgoType t = simplexMethodWithInitalization(m,n,nA,nb,nc,nx,val,iindb,indn);
	indb.clear();indb.reserve(iindb.size());
	for ( int i = 0 ; i < iindb.size() ; ++ i ){
		if(iindb[i]<n)
			indb.push_back(iindb[i]);
	}
	solveCompleteBlockingSystem(A,b,indb,x);
	val = x.dot(c);
	return t;
}

}// End of namespace COPT

#endif