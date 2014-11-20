// Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
// Copyright (C) MathU

#ifndef SIMPLEX_HPP
#define SIMPLEX_HPP
namespace COPT
{


/*		The Linear Programming(LP) problem is assuemd to be with
 *		standard form. The constraint contains two parts: Ax=b and x>=0.
 *		Any non-standard problem will be transformed at first.
 *		The very basic algorithm is implemented. We use a two phase 
 *		algorithm. The first phase is used to find a feasible x.
 */
template<class kernel>
class SimplexSolver
{
private:
	typedef typename kernel::ScalarType			ScalarType;
	typedef typename kernel::Vector				Vector;
	typedef typename kernel::Matrix				Matrix;
	// typedef LinearEqualConstraint<kernel>		ECons;


	/** the variable of the problem*/
	//%{
	/** the equal constraint */
	// ECons 				__econ;
	size_t 				__n;
	size_t 				__m;
	Vector				__x;
	Vector				__c;

	// type of the result of a step
	enum StepType{
		Normal,
		Optimal,
		Unbounded
	};

	enum AlgoType{
		Success,
		Unbound,
		NotFeasible
	};
	//%}

public:

	/** static functions of simplex solver */
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
						const std::vector<size_t>& indb,
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
						const std::vector<size_t>& indb,
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
 						std::vector<size_t>& indb,
 						std::vector<size_t>& indn);

 	static inline bool findPivoting(
 						const Vector& xb,
 						const Vector& d,
 						const int q,
 						std::vector<size_t>& indb,
 						std::vector<size_t>& indn);

 	/* 		The first phase of simplex method. The target is to find the feasible
	 *		point of the problem. 
	 *				min e^Tz
	 *				s.t. Ax+Ez=b
	 *					 (x,z)>=0
	 */
	static inline void findFeasiblePoint(
						const size_t m,
						const size_t n,
						const Matrix& A,
						const Vector& b,
						std::vector<size_t>& indb,
						std::vector<size_t>& indn,
						std::vector<size_t>& indz);
 	/* 		Simplex method algorithm when initial point is known.
 	 *		/param A:			equal constraint matrix
 	 *		/param b:			right hand vector of equal constraint
 	 *		/param c:			weight vector
 	 *		/param x:			initial feasible point on input and result on output
 	 */
 	static inline AlgoType simplexMethodWithInitalization(
 						const int m,
 						const int n,
 						const Matrix& A,
 						const Vector& b,
 						const Vector& c,
 						Vector& x,
 						std::vector<size_t>& indb,
 						std::vector<size_t>& indn);

 	static inline AlgoType simplexMethodWithIndices(
 						const int m,
 						const int n,
 						const Matrix& A,
 						const Vector& b,
 						const Vector& c,
 						Vector& x,
 						std::vector<size_t>& indb,
 						std::vector<size_t>& indn,
 						std::vector<size_t>& indz);
 	//%} end of static methods
}; // End of class SimplexSolver


template<class kernel>
void SimplexSolver<kernel>::solveBlockingSystem(
	const Matrix& A,
	const Vector& b,
	const std::vector<size_t>& indb,
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
	const std::vector<size_t>& indb,
	Vector& x)
{
	x.resize(A.cols());
	Vector xb;
	solveBlockingSystem(A,b,indb,xb);
	for ( int i = 0 ; i < indb.size() ; ++ i )
		x(indb[i]) = xb(i);
}


template<class kernel>
void SimplexSolver<kernel>::findFeasiblePoint(
	const size_t m,
	const size_t n,
	const Matrix& A,
	const Vector& b,
	std::vector<size_t>& indb,
	std::vector<size_t>& indn,
	std::vector<size_t>& indz)
{
	Matrix E = Matrix::identity(m,m),AE;
	for ( int i = 0 ; i < m ; ++ i ){
		if(b(i)<0)
			E(i,i) = static_cast<ScalarType>(-1.0);
		else
			E(i,i) = static_cast<ScalarType>(1.0);
	}
	Matrix::stCombineAlongColumn(A,E,AE);
	Vector e(m+n);
	for ( int i = 0 ; i < m ; ++ i )
		e(i+n)=1.0;
	Vector xz(m+n);
	for ( int i = 0 ; i < m ; ++ i )
		xz(i+n) = std::abs(b(i));
	std::vector<size_t> iindb(m),iindn(n); // the indices for B and N
	for ( size_t i = 0 ; i < m ; ++ i )
		iindb[i] = n+i;
	for ( size_t i = 0 ; i < n ; ++ i )
		iindn[i] = i;
	Vector x;
	simplexMethodWithInitalization(m,n,AE,b,e,x,iindb,iindn);

	// compute indb and indn. index less than n is kept
	indb.clear();indb.reserve(iindb.size());
	indn.clear();indn.reserve(iindn.size());
	indz.clear();indz.reserve(iindb.size());
	for ( size_t i = 0 ; i < iindb.size() ; ++ i )
	{
		if(iindb[i]<n)
			indb.push_back(iindb[i]);
		else
			indz.push_back(iindb[i]-n);
	}
	for (size_t i = 0 ; i < iindn.size() ; ++ i )
	{
		if(iindn[i]<n)
			indn.push_back(iindn[i]);
	}

}

template<class kernel>
typename SimplexSolver<kernel>::AlgoType SimplexSolver<kernel>::simplexMethodWithInitalization(
	const int m,
	const int n,
	const Matrix& A,
	const Vector& b,
	const Vector& c,
	Vector& x,
	std::vector<size_t>&indb,
	std::vector<size_t>&indn)
{
	StepType type;
	do{
		type = oneSimplexStep(A,b,c,indb,indn);
	}while(type == Normal);
	if ( type == Unbounded)
		return NotFeasible;
	else{
		Matrix B;
		B.columnBlockFromMatrix(A,indb);
		Vector xb = B.solve(b);
		x.resize(n);
		for ( size_t i = 0 ; i < m ; ++ i ){
			x(indb[i]) = xb(i);
		}
		solveCompleteBlockingSystem(A,b,indb,x);
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
	std::vector<size_t>&indb,
	std::vector<size_t>&indn,
	std::vector<size_t>&indz)
{
	// the difference is that some artificial variables still remains
	// the size is m - indb.size()
	size_t zsize = indz.size();
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
	std::vector<size_t> iindb(indb);
	for ( int i = 0 ; i < zsize ; ++ i )
		iindb.push_back(n+i);
	AlgoType t = simplexMethodWithInitalization(m,n,nA,nb,nc,nx,iindb,indn);
	indb.clear();indb.reserve(iindb.size());
	for ( int i = 0 ; i < iindb.size() ; ++ i ){
		if(iindb[i]<n)
			indb.push_back(iindb[i]);
	}
	solveCompleteBlockingSystem(A,b,indb,x);
	return t;
}

/** implementation part*/
template<class kernel>
bool SimplexSolver<kernel>::findPivoting(const Vector& xb,const Vector& d,const int q,std::vector<size_t>& indb,std::vector<size_t>& indn)
{
	int p = -1;
	ScalarType min = INFTY;
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
	std::vector<size_t>& indb,
	std::vector<size_t>& indn
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

}// End of namespace COPT

#endif