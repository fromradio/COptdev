//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef MatrixBase_H
#define MatrixBase_H

/*
 *				class 'MatrixBase'
 *					dervied from base class 'Array'
 */
namespace COPT
{
/*
 *	Class of 'MatrixBase'
 *		the data is stored column by column
 *		the matrix is assumed not to be transpose or symmetric
 *		but symmetric and transpose flag is designed
 */
template<class FT,class I = int>
class MatrixBase
	:
	public Array<FT,I>
{
public:

	/** the scalar type */
	typedef 			FT 				 			scalar;
	/** pod scalar type */
	typedef typename get_pod_type<scalar>::type 	podscalar;
	/** the size type used */
	typedef 			I 							index;						
	/** define fthe category */
	typedef 			matrix_object				ObjectCategory;
	/** define the trait */
	typedef 			KernelTrait<FT,index>		Kernel;

private:
	/** definition used in implementation */
	typedef 			VectorBase<FT,index>			Vector;
	typedef 			Array<FT,index>					Arr;

	/**			private variables			*/
	//%{

	/** the size of rows */
	index					__rows;

	/** the size of columns */
	index 					__cols;

	/** whether the matrix is symmetric */
	bool					__sym;

	/** whether the matrix is transpose */
	bool 					__trans;

	//%}
public:
	/** constructor and deconstructor */
	//%{
	/** default constructor */
	MatrixBase();

	MatrixBase(const index m, const index n, const scalar* data=NULL);

	/** Copy assignment */
	MatrixBase(const MatrixBase& mat);

	/** deconstructor */
	~MatrixBase();
	//%} end of constructor and deconstructor



	/**	getters and setters*/
	//%{
	/** get the number of rows */
	const index&		rows() const;

	/** get the number of columns */
	const index& 		cols() const;

	/**	matlab-like element getter */
	scalar& operator() ( const index i , const index j);

	const scalar& operator() ( const index i, const index j) const;

	/** get the element using Arr */
	const scalar& data( const index i ) const;

	/**	obtain the i-th row */
	Vector row( const index num );
	const Vector row( const index num ) const;
	/**	obtain the i-th column */
	Vector col( const index num );
	const Vector col( const index num ) const;


	/** set element using Arr */

	void set ( const index i , const scalar value );


	/** resize the matrix */
	void resize ( index m , index n );

	/** set the matrix to be symmetric */
	void setSymmetricFlag( bool sym );
	bool isSymmetric() const;

	/** set the matrix to be transpose mode */
	void setTransposeFlag( bool sym );
	bool isTranspose() const;

	/** Copy operation */
	MatrixBase& operator= ( const MatrixBase& mat );

	/*
		Mathematical operations
	*/

	/** summation */
	MatrixBase operator+ (const MatrixBase& mat);

	// subtraction
	// need to be tested
	MatrixBase operator- (const MatrixBase& mat);

	/** matrix multiplications */
	VectorBase<scalar,index> operator* ( const VectorBase<scalar,index>& vec ) const;

	/** matrix and matrix multiplication */
	MatrixBase operator* ( const MatrixBase& mat ) const;


	// multiplication between a scalar and a matrix
	friend MatrixBase operator* (const scalar s,const MatrixBase& mat)
	{
		MatrixBase result(mat.rows(),mat.cols());
		for ( index i = 0 ; i < mat.rows() ; ++ i )
			for ( index j = 0 ; j < mat.cols() ; ++ j )
				result(i,j) = mat(i,j)*s;
		return result;
	}

	// multiplication between a scalar and a matrix
	friend MatrixBase operator* (const MatrixBase& mat,const scalar s)
	{
		return s*mat;
	}

	// transpose
	MatrixBase transpose() const;

	/** transpose muliplication */
	Vector transMulti( const Vector& vec ) const;

	MatrixBase transMulti(const MatrixBase& mat) const;

	/*				overload of ostream
	 */
	friend std::ostream& operator<<(std::ostream& os,const MatrixBase& mat)
	{
		for ( index i = 0 ; i < mat.rows() ; ++ i ){
			for ( index j = 0 ; j < mat.cols() ; ++ j )
				os<<mat(i,j)<<' ';
			os<<std::endl;
		}
		return os;
	}

	/*				solve linear system
	 *
	 */
	 Vector lapackSolve( const Vector& vec ){
		if (__rows != __cols )
			throw COException("Solving Error: the matrix is not square!");
		if (__rows != vec.size() )
			throw COException("Solving Error: the size of right hand vector is wrong!");
		Vector v(vec);
		int pivot[__rows], c2 = 1,info;
		scalar* a = new scalar[__rows*__cols];
		blas::copt_blas_copy(__rows*__cols,this->dataPtr(),1,a,1);
		copt_lapack_gesv(&__rows,&c2,a,
			&__rows,pivot,v.dataPtr(),
			&__rows,
			&info);
		delete[] a; // delete the temporaray array
		return v;
	}
	
// #ifdef EIGEN
// 	VectorBase<scalar,index> solve(const VectorBase<scalar,index>& vec){
// 		// currently we use eigen to solve it
// 		Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic> matrix(__rows,__cols);
// 		for ( index i = 0 ; i < __rows ; ++ i )
// 		{
// 			for ( index j = 0 ;  j < __cols ; ++ j )
// 			{
// 				matrix(i,j) = this->operator()(i,j);
// 			}
// 		}
// 		Eigen::Matrix<scalar,Eigen::Dynamic,1> vector(vec.size());
// 		for ( index i = 0 ; i < vec.size() ; ++ i )
// 		{
// 			vector(i) = vec[i];
// 		}
// 		Eigen::Matrix<scalar,Eigen::Dynamic,1> result = matrix.colPivHouseholderQr().solve(vector);
// 		return VectorBase<scalar,index>(result);
// 	}
	Vector solve( const Vector& vec ){
		return lapackSolve(vec);
	}

// #ifdef USE_LAPACK
	Vector leastSquareSolve( const Vector& vec){
		if( __cols != vec.size() )
			throw COException("least square error: the size is not consistent!");
		scalar *a = new scalar[__rows*__cols];
		scalar *b = new scalar[vec.size()];
		scalar *s = new scalar[std::min(__rows,__cols)];
		index info = -1;
		index rank;
		blas::copt_blas_copy(__rows*__cols,this->dataPtr(),1,a,1);
		blas::copt_blas_copy(vec.size(),vec.dataPtr(),1,b,1);
		copt_lapack_gelss(__rows,__cols,1,a,__rows,b,vec.size(),s,-1.0,&rank,&info);
		Vector result(__cols,b);
		delete[]a;
		delete[]b;
		return result;
	}
// #endif
	


	/*			Special MatrixBase
	 *
	 */
	static MatrixBase identity(index m,index n){
		MatrixBase result(m,n);
		// std::cout<<result<<std::endl;
		index min = std::min(m,n);
		for ( index i = 0 ; i < min ; ++ i )
			result(i,i) = static_cast<scalar>(1.0);
		return result;
	}

	static MatrixBase identity(index m,index n,const scalar s);

	/** Blocking methods */
	//%{
	
	/** Blocking matrix from a given matrix with specific row numbers and column numbers*/
	void blockFromMatrix(
		const MatrixBase& mat,
		const std::set<index>& rownums,
		const std::set<index>& colnums);
	
	/** Blocking matrix from a given matrix with just columns */
	void columnBlockFromMatrix(
		const MatrixBase& mat,
		const std::set<index>& colnums);
	
	/** Blocking matrix from a given matrix with just rows */
	void rowBlockFromMatrix(
		const MatrixBase& mat,
		const std::set<index>& rownums);

	/** Blocking matrix from a given matrix, order is not considered */
	void blockFromMatrix(
		const MatrixBase& mat,
		const std::vector<index>& rownums,
		const std::vector<index>& colnums);

	/** blocking matrix using iterator */
	template<class InputIterator>
	void blockFromMatrix(
		const MatrixBase& mat,
		const InputIterator& rowbegin,
		const InputIterator& rowend,
		const InputIterator& colbegin,
		const InputIterator& colend );

	/** Blocking matrix from a given matrix, order is not considered */
	void columnBlockFromMatrix(
		const MatrixBase& mat,
		const std::vector<index>& colnums);

	/** column blocking from a given matrix */
	template<class InputIterator>
	void columnBlockFromMatrix(
		const MatrixBase& mat,
		const InputIterator& colbegin,
		const InputIterator& colend);

	/** Blocking matrix from a given matrix, order is not considered */
	void rowBlockFromMatrix(
		const MatrixBase& mat,
		const std::vector<index>& rownums);

	template<class InputIterator>
	void rowBlockFromMatrix(
		const MatrixBase& mat,
		const InputIterator& rowbegin,
		const InputIterator& rowend);

	//%}

	/** combining methods */
	//%{

	/** combining along row direction */
	void combineAlongRow(
		const MatrixBase& m1,
		const MatrixBase& m2);
	/** combining along column direction */
	void combineAlongColumn(
		const MatrixBase& m1,
		const MatrixBase& m2);

	/** combine along row direction taking matrix as parameter */
	static inline void stCombineAlongRow(
		const MatrixBase& m1,
		const MatrixBase& m2,
		MatrixBase& m);
	/** combine along column direction taking matrix as parameter */
	static inline void stCombineAlongColumn(
		const MatrixBase&m1,
		const MatrixBase& m2,
		MatrixBase& m);
	//%}

	/** set a random matrix */
	void setRandom( const index rows, const index cols);
	static inline MatrixBase random( const index rows, const index cols);
	/** compute A^TA of a given matrix */
	void mtm(MatrixBase& mat) const;

	/** norms */
	//%{
	/** compute the operation norm */
	podscalar operationNorm() const;
	//%}

};// End of class MatrixBase


/*
 *		class Triplet for sparse matrix assignment
 */
template<class Scalar,class I>
struct TripletBase
{

public:
	typedef Scalar 								scalar;
	typedef I 									index;
	typedef KernelTrait<Scalar,I>				kernel;
	typedef typename kernel::size 				size;
private:
	/** private variables */
	//%{

	/** row index */
	index				__r;

	/** column index */
	index 				__c;

	/** value */
	Scalar 				__v;

	//%}

	/** private constructors */
	//%{
	TripletBase();
	//%}

public:

	/** constructor and deconstructor */
	//%{

	/** the only constructor */
	TripletBase(
		const index r,
		const index c,
		const Scalar v);

	~TripletBase();
	//%}

	/** getter (no other setter is allowed) */
	//%{

	/** row index */
	const index& rowIndex() const;

	/** column index */
	const index& columnIndex() const;

	/** value */
	const scalar& value() const;
	//%}
};

/*		compare two triplets according to row index
 *
 */
template<class Triplet>
struct rowComparison
{
	bool operator()(const Triplet& t1,const Triplet& t2);
};

/*		compare two triplets accordSize to column index
 *
 */
template<class Triplet>
struct columnComparison
{
	bool operator()(const Triplet& t1,const Triplet& t2);
};

/*		Sparse matrix class
 *		the sparse matrix is designed for solve sparse linear systems
 */

template<class SpMatrix>
class UMFLinearSolver;

template<class T,class I>
class SpMatrixBase
{
public:
	typedef 	T 						scalar;
	typedef  	I 	 					index;
	typedef 	TripletBase<T,I>		Triplet;
	typedef 	sp_matrix_object 		ObjectCategory;
	typedef 	KernelTrait<T,I>		kernel;
private:

	typedef 	VectorBase<T,I>			Vector;
	/** private variables */
	//%{

	/** the number of rows */
	index 				__rows;

	/** the number of columns */
	index 				__cols;

	/** the number of elements */
	index 				__elesize;

	/** the col poSizeers */
	index*				__colptr;

	/** the indices of the rows */
	index*				__rowind;

	/** the values */
	scalar*		 		__vals;

	/** static zero */
	static const scalar __zero;

	//%}

	/** private functions */
	//%{

	/** judge the rationality */
	void judgeRationality();

	//%}

public:

	/** constructor and deconstructor */
	//%{

	/** default constructor */
	SpMatrixBase();

	SpMatrixBase(
		const index 				rows,
		const index 				cols,
		const index 				elesize,
		const index*			 	colptr,
		const index*				rowind,
		const scalar*				vals);

	SpMatrixBase(
		const SpMatrixBase& );

	/** deconstructor */
	~SpMatrixBase();
	//*}

	/** getter and setter */
	//%{

	/** traditional setter of sparse matrix*/
	void setSparseMatrix(
		const index 					rows,
		const index 					cols,
		const index 					elesize,
		const index*			 		colptr,
		const index*			 		rowind,
		const scalar*			 		vals);

	/** overload of operator = */
	SpMatrixBase& operator = (const SpMatrixBase& );

	/** set from triplets */
	void setFromTriplets(
		const index rows,
		const index cols,
		std::vector<Triplet>& triplets);

	/** set from triplet iterator */
	template<class InputIterator>
	void setFromTriplets(
		const index rows,
		const index cols,
		const InputIterator& begin,
		const InputIterator& end);

	/** fast setting from triplets iterator (like mtx file) */
	template<class InputIterator>
	void fastSetFromTriplets(
		const index rows,
		const index cols,
		const InputIterator& begin,
		const InputIterator& end);

	/** clear the data */
	void clear();

	/** get the row number */
	const index& rows() const;

	/** get the column number */
	const index& cols() const;

	/** get the element size */
	const index& elementSize() const;

	/** get the column poSizeer */
	const index* columnPointer() const; 

	/** get the row indices */
	const index* rowIndex() const;

	/** get the values */
	const scalar* values() const;

	/** scale with a scalar s */
	void scale(const scalar s);

	/** negative sign of the matrix */
	void neg();

	//%}

	/** element access */
	//%{

	/** only const access is allowed */
	const scalar& operator()(
		const index i,
		const index j) const;

	//%}

	/** operations */
	//%{

	/** summation */
	SpMatrixBase operator+(const SpMatrixBase& mat) const;

	/** subtraction */
	SpMatrixBase operator-(const SpMatrixBase& mat) const;

	/** negative sign */
	SpMatrixBase operator-() const;

	/** multiplication with another sparse matrix */
	SpMatrixBase operator*(const SpMatrixBase& mat) const;

	/** multiplication with vector */
	Vector operator*(const Vector& vec) const;

	/** multiplication with scalar */
	SpMatrixBase operator*(const scalar s) const;

	/** transform to a dense matrix */
	MatrixBase<scalar,index> toDenseMatrix() const;

	/** solve a linear system */
	VectorBase<scalar,index> solve(const VectorBase<scalar,index>& vec);

	//%}

}; // end of class SpMatrixBase


/** Sparse matrix related operator */
template<class scalar,class index>
SpMatrixBase<scalar,index> operator* (const scalar s,const SpMatrixBase<scalar,index>& mat);
template<class scalar,class index,class T>
SpMatrixBase<scalar,index> operator* (const T s,const SpMatrixBase<scalar,index>& mat);


}// End of namespace COPT
#endif