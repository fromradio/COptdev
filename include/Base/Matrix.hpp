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
	Class of 'MatrixBase'
		the data is stored column by column
*/
template<class FT,class S = longsize>
class MatrixBase
	:
	public Array<FT,S>
{
public:

	/** the scalar type */
	typedef 			FT 				 			ScalarType;
	/** the size type used */
	typedef 			S 							Size;						
	/** define fthe category */
	typedef 			matrix_tag 					Category;
	/** define the trait */
	typedef 			KernelTrait<FT,Size>		Kernel;

private:
	/** definition used in implementation */
	typedef 			VectorBase<FT,Size>			Vector;
	typedef 			Array<FT,Size>				Arr;

	/**			private variables			*/
	//%{

	/** the size of rows */
	Size					__rows;

	/** the size of columns */
	Size 					__cols;

	//%}
public:
	/** constructor and deconstructor */
	//%{
	/** default constructor */
	MatrixBase();

	MatrixBase(const Size m, const Size n,ScalarType* data=NULL);

	/** Copy assignment */
	MatrixBase(const MatrixBase& mat);

	/** deconstructor */
	~MatrixBase();
	//%} end of constructor and deconstructor



	/**	getters and setters*/
	//%{
	/** get the number of rows */
	const Size&		rows() const;

	/** get the number of columns */
	const Size& 		cols() const;

	/**	matlab-like element getter */
	ScalarType& operator() ( const Size i , const Size j);

	const ScalarType& operator() ( const Size i, const Size j) const;

	/** get the element using Arr */
	const ScalarType& data( const Size i ) const;

	/**	obtain the i-th row */
	Vector row( const Size num );
	const Vector row( const Size num ) const;
	/**	obtain the i-th column */
	Vector col( const Size num );
	const Vector col( const Size num ) const;


	// set element using Arr

	void set ( const Size i , ScalarType value ){
		if ( i < 0 )
			throw COException("MatrixBase error: index is less that zero!");
		else if ( i >= __rows*__cols )
			throw COException("MatrixBase error: index is out of range!");
		else
			this->operator[](i) = value;
	}

	/** resize the matrix */
	void resize ( Size m , Size n );

	/*
		Copy operation
	*/
	MatrixBase& operator= ( const MatrixBase& mat ) {
		if( __rows != mat.rows() || __cols != mat.cols() ){
			__rows = mat.rows();
			__cols = mat.cols();
			SAFE_DELETE_ARRAY(this->dataPtr());
			this->reset(__rows*__cols);
		}

		for ( Size i = 0 ; i < __rows*__cols ; ++ i )
			this->operator[](i) = mat.data(i);
		return *this;
	}

	/*
		Mathematical operations
	*/

	// summation
	// need to be tested
		
	MatrixBase operator+ (const MatrixBase& mat) {
		if ( __rows != mat.rows() || __cols != mat.cols() ) 
			throw COException("MatrixBase summation error: the size of two matrices are not consistent!");
		MatrixBase result(__rows,__cols);
		for ( Size i = 0 ; i < __rows*__cols ; ++ i )
			result.set(i,this->operator[](i)+mat.data(i));
		return result;
	}

	// subtraction
	// need to be tested
	MatrixBase operator- (const MatrixBase& mat) {
		if ( __rows != mat.rows() || __cols != mat.cols() ) 
			throw COException("MatrixBase subtraction error: the size of two matrices are not consistent!");
		MatrixBase result(__rows,__cols);
		for ( Size i = 0 ; i < __rows*__cols ; ++ i )
			result.set(i,this->operator[](i)-mat.data(i));
		return result;
	}

	// multiply
	// need to be tested
	VectorBase<ScalarType,Size> operator* ( const VectorBase<ScalarType,Size>& vec ) const{
		if ( __cols != vec.size() )
			throw COException("MatrixBase multiply error: the size of MatrixBase and vector are not consistent!");
		VectorBase<ScalarType,Size> result(__rows);
		for ( Size i = 0 ; i < __rows ; ++ i ){
			for ( Size j = 0 ; j < __cols ; ++ j )
				result[i]+= operator()(i,j)*vec[j];
		}
		return result;
	}
	// need to be tested
	MatrixBase operator* ( const MatrixBase& mat ) const{
		if ( __cols != mat.rows() )
			throw COException("MatrixBase multiply error: the size of two matrices are not consistent!");
		MatrixBase result (__rows,mat.cols());
		for ( Size i = 0 ; i < __rows ; ++ i )
			for ( Size j = 0 ; j < mat.cols() ; ++ j )
				for ( Size k = 0 ; k < __cols ; ++ k )
					result(i,j) += operator()(i,k)*mat(k,j);
		return result;
	}


	// multiplication between a scalar and a matrix
	friend MatrixBase operator* (const ScalarType s,const MatrixBase& mat)
	{
		MatrixBase result(mat.rows(),mat.cols());
		for ( Size i = 0 ; i < mat.rows() ; ++ i )
			for ( Size j = 0 ; j < mat.cols() ; ++ j )
				result(i,j) = mat(i,j)*s;
		return result;
	}

	// multiplication between a scalar and a matrix
	friend MatrixBase operator* (const MatrixBase& mat,const ScalarType s)
	{
		return s*mat;
	}

	// transpose
	MatrixBase transpose() const{
		MatrixBase result(__cols,__rows);
		for ( Size i = 0 ; i < __cols ; ++ i )
			for ( Size j = 0 ; j < __rows ; ++ j )
				result(i,j) = this->operator()(j,i);
		return result;
	}

	/*				overload of ostream
	 */
	friend std::ostream& operator<<(std::ostream& os,const MatrixBase& mat)
	{
		for ( Size i = 0 ; i < mat.rows() ; ++ i ){
			for ( Size j = 0 ; j < mat.cols() ; ++ j )
				os<<mat(i,j)<<' ';
			os<<std::endl;
		}
		return os;
	}

	/*				solve linear system
	 *
	 */
#ifdef EIGEN
	VectorBase<ScalarType,Size> solve(const VectorBase<ScalarType,Size>& vec){
		// currently we use eigen to solve it
		Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic> matrix(__rows,__cols);
		for ( Size i = 0 ; i < __rows ; ++ i )
		{
			for ( Size j = 0 ;  j < __cols ; ++ j )
			{
				matrix(i,j) = this->operator()(i,j);
			}
		}
		Eigen::Matrix<ScalarType,Eigen::Dynamic,1> vector(vec.size());
		for ( Size i = 0 ; i < vec.size() ; ++ i )
		{
			vector(i) = vec[i];
		}
		Eigen::Matrix<ScalarType,Eigen::Dynamic,1> result = matrix.colPivHouseholderQr().solve(vector);
		return VectorBase<ScalarType,Size>(result);
	}
#endif


	/*			Special MatrixBase
	 *
	 */
	static MatrixBase identity(Size m,Size n){
		MatrixBase result(m,n);
		// std::cout<<result<<std::endl;
		Size min = std::min(m,n);
		for ( Size i = 0 ; i < min ; ++ i )
			result(i,i) = static_cast<ScalarType>(1.0);
		return result;
	}

	static MatrixBase identity(Size m,Size n,const ScalarType s);

	/** Blocking methods */
	//%{
	
	/** Blocking matrix from a given matrix with specific row numbers and column numbers*/
	void blockFromMatrix(
		const MatrixBase& mat,
		const std::set<Size>& rownums,
		const std::set<Size>& colnums);
	
	/** Blocking matrix from a given matrix with just columns */
	void columnBlockFromMatrix(
		const MatrixBase& mat,
		const std::set<Size>& colnums);
	
	/** Blocking matrix from a given matrix with just rows */
	void rowBlockFromMatrix(
		const MatrixBase& mat,
		const std::set<Size>& rownums);

	/** Blocking matrix from a given matrix, order is not considered */
	void blockFromMatrix(
		const MatrixBase& mat,
		const std::vector<Size>& rownums,
		const std::vector<Size>& colnums);

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
		const std::vector<Size>& colnums);

	/** column blocking from a given matrix */
	template<class InputIterator>
	void columnBlockFromMatrix(
		const MatrixBase& mat,
		const InputIterator& colbegin,
		const InputIterator& colend);

	/** Blocking matrix from a given matrix, order is not considered */
	void rowBlockFromMatrix(
		const MatrixBase& mat,
		const std::vector<Size>& rownums);

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

};// End of class MatrixBase


/*
 *		class Triplet for sparse matrix assignment
 */
template<class Scalar,class Size = size_t>
struct TripletBase
{

public:
	typedef Scalar 				ScalarType;
private:
	/** private variables */
	//%{

	/** row index */
	Size				__r;

	/** column index */
	Size 				__c;

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
		const Size r,
		const Size c,
		const Scalar v);

	~TripletBase();
	//%}

	/** getter (no other setter is allowed) */
	//%{

	/** row index */
	const Size& rowIndex() const;

	/** column index */
	const Size& columnIndex() const;

	/** value */
	const ScalarType& value() const;
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

template<class T,class S = longsize>
class SpMatrixBase
{
public:
	typedef 	T 						ScalarType;
	typedef  	S 	 					Size;
	typedef 	TripletBase<T,S>		Triplet;
	typedef 	matrix_tag 				Category;
private:

	typedef 	VectorBase<T,S>			Vector;
	/** private variables */
	//%{

	/** the number of rows */
	Size 				__rows;

	/** the number of columns */
	Size 				__cols;

	/** the number of elements */
	Size 				__elesize;

	/** the col poSizeers */
	Size*				__colptr;

	/** the indices of the rows */
	Size*				__rowind;

	/** the values */
	ScalarType*		 	__vals;

	/** static zero */
	static const ScalarType __zero;

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
		const Size 				rows,
		const Size 				cols,
		const Size 				size,
		const Size*			 	colptr,
		const Size*				rowind,
		const ScalarType*			vals);

	SpMatrixBase(
		const SpMatrixBase& );

	/** deconstructor */
	~SpMatrixBase();
	//*}

	/** getter and setter */
	//%{

	/** traditional setter of sparse matrix*/
	void setSparseMatrix(
		const Size 					rows,
		const Size 					cols,
		const Size 					size,
		const Size*			 		colptr,
		const Size*			 		rowind,
		const ScalarType*			 	vals);

	/** overload of operator = */
	SpMatrixBase& operator = (const SpMatrixBase& );

	/** set from triplets */
	void setFromTriplets(
		const Size rows,
		const Size cols,
		std::vector<Triplet>& triplets);

	/** set from triplet iterator */
	template<class InputIterator>
	void setFromTriplets(
		const Size rows,
		const Size cols,
		const InputIterator& begin,
		const InputIterator& end);

	/** clear the data */
	void clear();

	/** get the row number */
	const Size& rows() const;

	/** get the column number */
	const Size& cols() const;

	/** get the element size */
	const Size& elementSize() const;

	/** get the column poSizeer */
	const Size* columnPointer() const; 

	/** get the row indices */
	const Size* rowIndex() const;

	/** get the values */
	const ScalarType* values() const;

	/** scale with a scalar s */
	void scale(const ScalarType s);

	/** negative sign of the matrix */
	void neg();

	//%}





	/** element access */
	//%{

	/** only const access is allowed */
	const ScalarType& operator()(
		const Size i,
		const Size j) const;

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
	SpMatrixBase operator*(const ScalarType s) const;

	/** transform to a dense matrix */
	MatrixBase<ScalarType,Size> toDenseMatrix() const;

	/** solve a linear system */
	VectorBase<ScalarType,Size> solve(const VectorBase<ScalarType,Size>& vec);

	//%}

}; // end of class SpMatrixBase


/** Sparse matrix related operator */
template<class ScalarType,class Size>
SpMatrixBase<ScalarType,Size> operator* (const ScalarType s,const SpMatrixBase<ScalarType>& mat);
template<class ScalarType,class Size,class T>
SpMatrixBase<ScalarType,Size> operator* (const T s,const SpMatrixBase<ScalarType>& mat);


}// End of namespace COPT
#endif