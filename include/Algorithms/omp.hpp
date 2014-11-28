//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef OMP_H
#define OMP_H

namespace COPT{
/*			This file describes a famous algorithm called orthgonal
 *			matching pursuit (OMP). The algorithm solves sparse problem:
 *				min \|Ax-b\|_2 s.t.\|x\|_0 \leq s
 *			Input:
 *				Atom or dictionary matrix A
 *				Observation x
 *				sparsity s
 *			Output:
 *				the coefficient vector
 */
template<class kernel>
class OMPSolver
	:
	noncopyable
{
private:
	typedef typename kernel::ScalarType 	Scalar;
	typedef typename kernel::Size 			Size;
	typedef typename kernel::Vector 		Vector;
	typedef typename kernel::Matrix 		Matrix;

	/** 		private variables */
	//%{

	/** the problem is with size m*n */
	Size 					__m;
	Size 					__n;

	/** atom or dictionary matrix A */
	Matrix 					__A;

	/** the transpose of A */
	Matrix 					__AT;

	/** A^T*A */
	Matrix 					__ATA;

	/** obeservation vector b */
	Vector 			 		__b;

	/** sparsity */
	Size 					__s;

	/** the result */
	Vector 					__x;

	/** the indices of using atoms */
	std::list<Size> 		__used_indices;

	/** the vector of residual */
	Vector 					__residual;

	/** the coefficient vector */
	Vector 					__coeff;

	/** whether the atom matrix is normalized*/
	bool 					__is_normalized;

	/** storing the norm of very column is input A is not normalized */
	Vector 					__norms;

	/** whether ATA exists */
	bool 					__is_batch;

	/** threshold */
	Scalar 					__error;
	//%}

	OMPSolver();

	/** 	whether zero column exists which is not allowed in the algorithm */
	inline bool		checkZeroColumn();

	/** 	normalize the column of argument */
	static inline void 	normalizeAtomMatrix( 
		Matrix& A , 
		Vector& norms );

	/**		normalize the column of the argument and the big D */
	static inline void 	normalizeAtomMatrix( 
		Matrix& A , 
		Matrix& ATA , 
		Vector& norms );

	/**		find the next index */
	static inline Size 	findIndex( 
		const Matrix& A, 
		const Vector& residual );

	/** 	update the coefficient of the OMP solver */
	static inline void 	updateCoefficient( 
		const Matrix& A,
		const Matrix& AT,
		const Vector& b, 
		const std::list<Size>& indices, 
		Vector& coeff ,
		Vector& residual);

	/*		Update the coefficient of OMP solver.
	 *		Batch OMP algorithm is used.
	 */
	static inline void 	updateCoefficient( 
		const Matrix& A,
		const Matrix& AT,
		const Matrix& ATA,
		const Vector& b, 
		const std::list<Size>& indices, 
		Vector& coeff ,
		Vector& residual);

	/**		kernel iterations of OMP algorithm */
	inline void 		ompSolve( 
		const Vector& b);

public:

	/**			constructors and deconstructors 		*/
	//%{

	OMPSolver(const Matrix&A,const bool isnormalized = false);

	OMPSolver(const Matrix&A,const Matrix& ATA, const bool isnormalized = false );

	//%}



	/**			    getter and setter 					*/
	//%{

	/** set the sparsity */
	void setSparsity(const Size s);

	/** get the sparsity */
	const Size& sparsity() const;

	/** get the matrix */
	const Matrix& atomMatrix() const;

	/** get the observation */
	const Vector& observation() const;

	/** get the result */
	const Vector& result() const;

	/** get the final residual */
	const Vector& finalResidual() const;

	/** get the fitting error */
	const Scalar& fittingError() const;

	//%}

	/**				solver interface 					*/
	void solve(const Vector& b);
	void solve(const Vector& b,const Size s);

};

/** 		Implementation of OMPSolver				*/

template<class kernel>
bool OMPSolver<kernel>::checkZeroColumn( )
{
	for ( Size i = 0 ; i < __n ; ++ i )
	{
		if(IS_ZERO(__A.col(i).squaredNorm()))
			return true;
	}
	return false;
}

template<class kernel>
void OMPSolver<kernel>::normalizeAtomMatrix( Matrix& A,Vector& norms)
{
	norms.resize(A.cols());
	for ( Size i = 0 ; i < A.cols() ; ++ i )
	{
		norms[i] = A.col(i).normalize();
	}
}

template<class kernel>
void OMPSolver<kernel>::normalizeAtomMatrix( Matrix& A, Matrix& ATA , Vector& norms )
{
	norms.resize(A.cols());
	for ( Size i = 0 ; i < A.cols() ; ++ i )
	{
		norms[i] = A.col(i).normalize();
	}
	for ( Size i = 0 ; i < A.cols() ; ++ i )
	{
		for ( Size j = 0 ; j < A.cols() ; ++ j )
		{
			ATA(i,j) = ATA(i,j)*(1.0/norms[i])*(1.0/norms[j]);
		}
	}
}

template<class kernel>
typename kernel::Size OMPSolver<kernel>::findIndex( const Matrix& A, const Vector& residual )
{
	// find the new column for orthogonal matching pursuit (OMP)
	// make sure that every column is normalized at first
	Size index = -1;
	Scalar maximal = -1.0;
	for ( Size i = 0 ; i < A.cols() ; ++ i )
	{
		Scalar v = std::abs(A.col(i).dot(residual));
		if ( v > maximal ){
			index = i;
			maximal = v;
		}
	}
	return index;
}

template<class kernel>
void OMPSolver<kernel>::updateCoefficient(const Matrix&A , const Matrix&AT, const Vector& b , const std::list<Size>& indices, Vector& coeff , Vector& residual )
{
	//		 not batch omp
	Matrix Ab,ATb;
	Ab.columnBlockFromMatrix(A,indices.begin(),indices.end());
	ATb.rowBlockFromMatrix(AT,indices.begin(),indices.end());
	coeff = (ATb*Ab).solve(ATb*b);
	residual = b - Ab*coeff;
}

template<class kernel>
void OMPSolver<kernel>::updateCoefficient(const Matrix&A , const Matrix&AT, const Matrix& ATA , const Vector& b ,const std::list<Size>& indices, Vector& coeff , Vector& residual )
{
	// 		batch omp
	Matrix Ab,ATb,Db;
	Ab.columnBlockFromMatrix(A,indices.begin(),indices.end());
	ATb.rowBlockFromMatrix(AT,indices.begin(),indices.end());
	Db.blockFromMatrix(ATA,indices.begin(),indices.end(),indices.begin(),indices.end());
	coeff = Db.solve(ATb*b);
	residual = b - Ab*coeff;
}

template<class kernel>
void OMPSolver<kernel>::ompSolve(const Vector& b)
{
	__used_indices.clear();
	__residual = b;
	if (!__is_batch)
	{
		for ( Size i = 0 ; i < __s ; ++ i )
		{
			__used_indices.push_back(findIndex(__A,__residual));
			updateCoefficient(__A,__AT,b,__used_indices,__coeff,__residual);
			if ( IS_ZERO(__residual.squaredNorm()) )
				break;
		}
	}
	else
	{
		for ( Size i = 0 ; i < __s ; ++ i )
		{
			__used_indices.push_back(findIndex(__A,__residual));
			updateCoefficient(__A,__AT,__ATA,b,__used_indices,__coeff,__residual);
			if ( IS_ZERO(__residual.squaredNorm()) )
				break;
		}
	}
	// set __x
	__x.resize(__n);
	Size i = 0;
	for (typename std::list<Size>::iterator iter = __used_indices.begin() ; iter != __used_indices.end() ; ++ iter )
	{
		__x(*iter) = __coeff(i++)*(1.0/__norms[*iter]);
	}
	__error = std::sqrt(__residual.squaredNorm());
}

template<class kernel>
OMPSolver<kernel>::OMPSolver(const Matrix& A,const bool isnormalizd)
	:
	__m(A.rows()),
	__n(A.cols()),
	__A(A),
	__AT(A.transpose()),
	__s(0),
	__is_normalized(isnormalizd),
	__is_batch(false)
{
	if(checkZeroColumn())
		throw COException("OMP Solver Error: there is zero column in atom matrix A. Solver ends!");
	if(!__is_normalized )
	{
		normalizeAtomMatrix(__A,__norms);
		__AT = __A.transpose();
	}
	else
	{
		__AT = __A.transpose();
	}
}

template<class kernel>
OMPSolver<kernel>::OMPSolver(const Matrix& A,const Matrix& ATA,const bool isnormalizd)
	:
	__m(A.rows()),
	__n(A.cols()),
	__A(A),
	__ATA(ATA),
	__s(0),
	__is_normalized(isnormalizd),
	__is_batch(true)
{
	if(checkZeroColumn())
		throw COException("OMP Solver Error: there is zero column in atom matrix A. Solver ends!");
	if(!__is_normalized ){
		normalizeAtomMatrix(__A,__ATA,__norms);
		__AT = __A.transpose();
	}
	else{
		__AT = __A.transpose();
	}
}

template<class kernel>
void OMPSolver<kernel>::setSparsity(const Size s)
{
	__s = s;
}

template<class kernel>
const typename kernel::Size& OMPSolver<kernel>::sparsity() const
{
	return __s;
}

template<class kernel>
const typename kernel::Matrix& OMPSolver<kernel>::atomMatrix() const
{
	return __A;
}

template<class kernel>
const typename kernel::Vector& OMPSolver<kernel>::observation() const
{
	return __b;
}

template<class kernel>
const typename kernel::Vector& OMPSolver<kernel>::result() const
{
	return __x;
}

template<class kernel>
const typename kernel::Vector& OMPSolver<kernel>::finalResidual() const
{
	return __residual;
}

template<class kernel>
const typename kernel::ScalarType& OMPSolver<kernel>::fittingError() const
{
	return __error;
}

template<class kernel>
void OMPSolver<kernel>::solve(const Vector& b){
	ompSolve(b);
}

template<class kernel>
void OMPSolver<kernel>::solve(const Vector& b,const Size s){
	setSparsity(s);
	ompSolve(b);
}

} // End of namespace COPT

#endif