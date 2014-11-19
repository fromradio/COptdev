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
	//%}

	/** private functions */
	//%{
	
 	
	/** the first phase */
	void phaseOne();
	//%}
public:
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

 	static inline bool findPivoting(const Vector& xb,const Vector& d,const int q,std::vector<size_t>& indb,std::vector<size_t>& indn);
}; // End of class SimplexSolver


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
	}
	if ( q == -1 )
		return Optimal;
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