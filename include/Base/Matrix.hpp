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
template<class FT>
class MatrixBase
	: public Array<FT>
{
public:

	// the scalar type
	typedef 			Array<FT>				Arr;
	typedef typename	Arr::ScalarType 		ScalarType;
	typedef 			VectorBase<FT>			Vector;


private:
	/** the size of rows */
	size_t					__rows;
	/** the size of columns */
	size_t 					__cols;
public:
	/** constructor and deconstructor */
	//%{
	/** default constructor */
	MatrixBase();

	MatrixBase(const size_t m, const size_t n,ScalarType* data=NULL);

	/** Copy assignment */
	MatrixBase(const MatrixBase& mat);

	/** deconstructor */
	~MatrixBase();
	//%} end of constructor and deconstructor



	/**	getters and setters*/
	//%{
	/** get the number of rows */
	const size_t&		rows() const;

	/** get the number of columns */
	const size_t& 		cols() const;

	/**	matlab-like element getter */
	ScalarType& operator() ( const int i , const int j);

	const ScalarType& operator() ( const int i, const int j) const;

	/** get the element using Arr */
	const ScalarType& data( const int i ) const;

	/**	obtain the i-th row */
	Vector row( const size_t num );
	/**	obtain the i-th column */
	Vector col( const size_t num );


	// set element using Arr

	void set ( int i , ScalarType value ){
		if ( i < 0 )
			throw COException("MatrixBase error: index is less that zero!");
		else if ( i >= __rows*__cols )
			throw COException("MatrixBase error: index is out of range!");
		else
			this->operator[](i) = value;
	}

	/** resize the matrix */
	void resize ( size_t m , size_t n );

	/*
		Copy operation
	*/
	MatrixBase& operator= ( const MatrixBase<ScalarType>& mat ) {
		if( __rows != mat.rows() || __cols != mat.cols() ){
			__rows = mat.rows();
			__cols = mat.cols();
			SAFE_DELETE_ARRAY(this->dataPtr());
			this->reset(__rows*__cols);
		}

		for ( int i = 0 ; i < __rows*__cols ; ++ i )
			this->operator[](i) = mat.data(i);
		return *this;
	}

	/*
		Mathematical operations
	*/

	// summation
	// need to be tested
		
	MatrixBase operator+ (const MatrixBase<ScalarType>& mat) {
		if ( __rows != mat.rows() || __cols != mat.cols() ) 
			throw COException("MatrixBase summation error: the size of two matrices are not consistent!");
		MatrixBase<ScalarType> result(__rows,__cols);
		for ( int i = 0 ; i < __rows*__cols ; ++ i )
			result.set(i,this->operator[](i)+mat.data(i));
		return result;
	}

	// subtraction
	// need to be tested
	MatrixBase operator- (const MatrixBase<ScalarType>& mat) {
		if ( __rows != mat.rows() || __cols != mat.cols() ) 
			throw COException("MatrixBase subtraction error: the size of two matrices are not consistent!");
		MatrixBase<ScalarType> result(__rows,__cols);
		for ( int i = 0 ; i < __rows*__cols ; ++ i )
			result.set(i,this->operator[](i)-mat.data(i));
		return result;
	}

	// multiply
	// need to be tested
	VectorBase<ScalarType> operator* ( const VectorBase<ScalarType>& vec ) const{
		if ( __cols != vec.size() )
			throw COException("MatrixBase multiply error: the size of MatrixBase and vector are not consistent!");
		VectorBase<ScalarType> result(__rows);
		for ( int i = 0 ; i < __rows ; ++ i ){
			for ( int j = 0 ; j < __cols ; ++ j )
				result[i]+= operator()(i,j)*vec[j];
		}
		return result;
	}
	// need to be tested
	MatrixBase<ScalarType> operator* ( const MatrixBase<ScalarType>& mat ) const{
		if ( __cols != mat.rows() )
			throw COException("MatrixBase multiply error: the size of two matrices are not consistent!");
		MatrixBase<ScalarType> result (__rows,mat.cols());
		for ( int i = 0 ; i < __rows ; ++ i )
			for ( int j = 0 ; j < mat.cols() ; ++ j )
				for ( int k = 0 ; k < __cols ; ++ k )
					result(i,j) += operator()(i,k)*mat(k,j);
		return result;
	}


	// multiplication between a scalar and a matrix
	friend MatrixBase<ScalarType> operator* (const ScalarType s,const MatrixBase<ScalarType>& mat)
	{
		MatrixBase<ScalarType> result(mat.rows(),mat.cols());
		for ( int i = 0 ; i < mat.rows() ; ++ i )
			for ( int j = 0 ; j < mat.cols() ; ++ j )
				result(i,j) = mat(i,j)*s;
		return result;
	}

	// multiplication between a scalar and a matrix
	friend MatrixBase<ScalarType> operator* (const MatrixBase<ScalarType>& mat,const ScalarType s)
	{
		return s*mat;
	}

	// transpose
	MatrixBase transpose() const{
		MatrixBase result(__cols,__rows);
		for ( int i = 0 ; i < __cols ; ++ i )
			for ( int j = 0 ; j < __rows ; ++ j )
				result(i,j) = this->operator()(j,i);
		return result;
	}

	/*				overload of ostream
	 */
	friend std::ostream& operator<<(std::ostream& os,const MatrixBase& mat)
	{
		for ( int i = 0 ; i < mat.rows() ; ++ i ){
			for ( int j = 0 ; j < mat.cols() ; ++ j )
				os<<mat(i,j)<<' ';
			os<<std::endl;
		}
		return os;
	}

	/*				solve linear system
	 *
	 */
#ifdef EIGEN
	VectorBase<ScalarType> solve(const VectorBase<ScalarType>& vec){
		// currently we use eigen to solve it
		Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic> matrix(__rows,__cols);
		for ( int i = 0 ; i < __rows ; ++ i )
		{
			for ( int j = 0 ;  j < __cols ; ++ j )
			{
				matrix(i,j) = this->operator()(i,j);
			}
		}
		Eigen::Matrix<ScalarType,Eigen::Dynamic,1> vector(vec.size());
		for ( int i = 0 ; i < vec.size() ; ++ i )
		{
			vector(i) = vec[i];
		}
		Eigen::Matrix<ScalarType,Eigen::Dynamic,1> result = matrix.colPivHouseholderQr().solve(vector);
		return VectorBase<ScalarType>(result);
	}
#endif


	/*			Special MatrixBase
	 *
	 */
	static MatrixBase identity(size_t m,size_t n){
		MatrixBase result(m,n);
		// std::cout<<result<<std::endl;
		size_t min = std::min(m,n);
		for ( int i = 0 ; i < min ; ++ i )
			result(i,i) = static_cast<ScalarType>(1.0);
		return result;
	}

	static MatrixBase identity(size_t m,size_t n,const ScalarType s);

	/** Blocking methods */
	//%{
	
	/** Blocking matrix from a given matrix with specific row numbers and column numbers*/
	void blockFromMatrix(
		const MatrixBase& mat,
		const std::set<size_t>& rownums,
		const std::set<size_t>& colnums);
	
	/** Blocking matrix from a given matrix with just columns */
	void columnBlockFromMatrix(
		const MatrixBase& mat,
		const std::set<size_t>& colnums);
	
	/** Blocking matrix from a given matrix with just rows */
	void rowBlockFromMatrix(
		const MatrixBase& mat,
		const std::set<size_t>& rownums);

	/** Blocking matrix from a given matrix, order is not considered */
	void blockFromMatrix(
		const MatrixBase& mat,
		const std::vector<size_t>& rownums,
		const std::vector<size_t>& colnums);

	/** Blocking matrix from a given matrix, order is not considered */
	void columnBlockFromMatrix(
		const MatrixBase& mat,
		const std::vector<size_t>& colnums);

	/** Blocking matrix from a given matrix, order is not considered */
	void rowBlockFromMatrix(
		const MatrixBase& mat,
		const std::vector<size_t>& rownums);

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
template<class Scalar>
struct TripletBase
{

public:
	typedef Scalar 				ScalarType;
private:
	/** private variables */
	//%{

	/** row index */
	size_t				__r;

	/** column index */
	size_t 				__c;

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
		const size_t r,
		const size_t c,
		const Scalar v);

	~TripletBase();
	//%}

	/** getter (no other setter is allowed) */
	//%{

	/** row index */
	const size_t& rowIndex() const;

	/** column index */
	const size_t& columnIndex() const;

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

/*		compare two triplets accordint to column index
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
template<class S>
class SpMatrixBase
{
public:
	typedef 	S 					ScalarType;
	typedef 	TripletBase<S>		Triplet;
private:

	typedef 	VectorBase<S>		Vector;
	/** private variables */
	//%{

	/** the number of rows */
	size_t 				__rows;

	/** the number of columns */
	size_t 				__cols;

	/** the number of elements */
	size_t 				__elesize;

	/** the indices of the rows */
	size_t*				__rowind;

	/** the col pointers */
	size_t*				__colptr;

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
		const size_t 				rows,
		const size_t 				cols,
		const size_t 				size,
		const size_t*				rowind,
		const size_t*			 	colptr,
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
		const size_t 					rows,
		const size_t 					cols,
		const size_t 					size,
		const size_t*			 		rowind,
		const size_t*			 		colptr,
		const ScalarType*			 	vals);

	/** overload of operator = */
	SpMatrixBase& operator = (const SpMatrixBase& );

	/** set from triplets */
	void setFromTriplets(
		const size_t rows,
		const size_t cols,
		std::vector<Triplet>& triplets);

	/** clear the data */
	void clear();

	/** get the row number */
	const size_t& rows() const;

	/** get the column number */
	const size_t& cols() const;

	/** get the element size */
	const size_t& elementSize() const;

	/** get the column pointer */
	const size_t* columnPointer() const; 

	/** get the row indices */
	const size_t* rowIndex() const;

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
		const size_t i,
		const size_t j) const;

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
	MatrixBase<ScalarType> toDenseMatrix() const;

	//%}

}; // end of class SpMatrixBase


/** Sparse matrix related operator */
template<class ScalarType>
SpMatrixBase<ScalarType> operator* (const ScalarType s,const SpMatrixBase<ScalarType>& mat);
template<class ScalarType,class T>
SpMatrixBase<ScalarType> operator* (const T s,const SpMatrixBase<ScalarType>& mat);


}// End of namespace COPT
#endif