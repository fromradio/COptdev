//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef MatrixBase_H
#define MatrixBase_H

// #include "Vector.h"
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
	typedef 			VectorBase<ScalarType>	Vector;


private:
	// the size of rows
	size_t					__rows;
	// the size of columns
	size_t 					__cols;
public:
	// default constructor
	MatrixBase()
		:
		Arr(),
		__rows(0),
		__cols(0)
	{

	}
	MatrixBase(size_t m,size_t n,ScalarType* data=NULL)
		:
		Arr(m*n,data),
		__rows(m),
		__cols(n)
	{
	}

	/*
		Copy assignment
	*/
	MatrixBase(const MatrixBase<ScalarType>& mat)
		:
		Arr(mat.rows()*mat.cols(),mat.dataPtr()),
		__rows(mat.rows()),
		__cols(mat.cols())
	{
	}

	~MatrixBase()
	{
	}

	/*
		basic getters
	*/
	// get the number of rows
	const size_t&		rows() const {return __rows;}
	// get the number of columns
	const size_t& 		cols() const {return __cols;}

	/*
		get the element of the MatrixBase
	*/

	ScalarType& operator() (int i , int j){
		if(i<0||j<0)
			throw COException("MatrixBase error: index is less than zero!");
		else if (i>=__rows||j>=__cols)
			throw COException("MatrixBase error: index is out of range!");
		else{
			return this->__data_ptr[j*__rows+i];
		}
	}
	const ScalarType& operator() (int i,int j) const {
		return const_cast<MatrixBase&>(*this).operator()(i,j);
	}

	// get the element using Arr

	const ScalarType& data( int i ) const{
		if ( i < 0 )
			throw COException("MatrixBase error: index is less that zero!");
		else if ( i >= __rows*__cols )
			throw COException("MatrixBase error: index is out of range!");
		else
			return this->__data_ptr[i];
	}

	/*
	 *			obtain the i-th row
	 */
	Vector row( int num ) const {
		if ( num < 0 )
			throw COException("MatrixBase error: row index is less that zero!");
		else if ( num >= __rows)
			throw COException("MatrixBase error: row index is out of range!");
		else{
			Vector result(__cols);
			for ( int i = 0 ; i < __cols ; ++ i ){
				result[i] = this->operator()(num,i);
			}
			return result;
		}
	}
	/*
	 *			obtain the i-th column
	 */
	Vector col( int num ) const {
		if ( num < 0 )
			throw COException("MatrixBase error: row index is less that zero!");
		else if ( num >= __cols)
			throw COException("MatrixBase error: row index is out of range!");
		else{
			Vector result(__rows);
			for ( int i = 0 ; i < __rows ; ++ i ){
				result[i] = this->operator()(i,num);
			}
			return result;
		}
	}


	// set element using Arr

	void set ( int i , ScalarType value ){
		if ( i < 0 )
			throw COException("MatrixBase error: index is less that zero!");
		else if ( i >= __rows*__cols )
			throw COException("MatrixBase error: index is out of range!");
		else
			this->__data_ptr[i] = value;
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
			SAFE_DELETE_ARRAY(this->__data_ptr);
			this->__data_ptr = new ScalarType[__rows*__cols];
		}

		for ( int i = 0 ; i < __rows*__cols ; ++ i )
			this->__data_ptr[i] = mat.data(i);
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
			result.set(i,this->__data_ptr[i]+mat.data(i));
		return result;
	}

	// subtraction
	// need to be tested
	MatrixBase operator- (const MatrixBase<ScalarType>& mat) {
		if ( __rows != mat.rows() || __cols != mat.cols() ) 
			throw COException("MatrixBase subtraction error: the size of two matrices are not consistent!");
		MatrixBase<ScalarType> result(__rows,__cols);
		for ( int i = 0 ; i < __rows*__cols ; ++ i )
			result.set(i,this->__data_ptr[i]-mat.data(i));
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
	void blockFromMatrix(const MatrixBase& mat,const std::set<size_t>& rownums,const std::set<size_t>& colnums);
	
	/** Blocking matrix from a given matrix with just columns */
	void columnBlockFromMatrix(const MatrixBase& mat,const std::set<size_t>& colnums);
	
	/** Blocking matrix from a given matrix with just rows */
	void rowBlockFromMatrix(const MatrixBase& mat,const std::set<size_t>& rownums);

	/** Blocking matrix from a given matrix, order is not considered */
	void blockFromMatrix(const MatrixBase& mat,const std::vector<size_t>& rownums,const std::vector<size_t>& colnums);
	/** Blocking matrix from a given matrix, order is not considered */
	void columnBlockFromMatrix(const MatrixBase& mat,const std::vector<size_t>& colnums);

	/** Blocking matrix from a given matrix, order is not considered */
	void rowBlockFromMatrix(const MatrixBase& mat,const std::vector<size_t>& rownums);
	//%}

	/** combining methods */
	//%{

	/** combining along row direction */
	void combineAlongRow(const MatrixBase& m1,const MatrixBase& m2);
	/** combining along column direction */
	void combineAlongColumn(const MatrixBase& m1,const MatrixBase& m2);

	/** combine along row direction taking matrix as parameter */
	static inline void stCombineAlongRow(const MatrixBase& m1,const MatrixBase& m2,MatrixBase& m);
	/** combine along column direction taking matrix as parameter */
	static inline void stCombineAlongColumn(const MatrixBase&m1,const MatrixBase& m2,MatrixBase& m);
	//%}

};// End of class MatrixBase

}// End of namespace COPT
#endif