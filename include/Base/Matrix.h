#ifndef MATRIX_H
#define MATRIX_H

// #include "Vector.h"
/*
 *				class 'Matrix'
 *					dervied from base class 'Array'
 */
namespace COPT
{
/*
	Class of 'Matrix'
		the data is stored column by column
*/
template<class FT>
class Matrix
	: public Array<FT>
{
public:

	// the scalar type
	typedef 			Array<FT>			Array;
	typedef typename	Array::ScalarType 	ScalarType;


private:
	// the size of rows
	int					__rows;
	// the size of columns
	int					__cols;
public:
	// default constructor
	Matrix()
		:
		Array(),
		__rows(0),
		__cols(0)
	{

	}
	Matrix(int m,int n,ScalarType* data=NULL)
		:
		Array(m*n,data),
		__rows(m),
		__cols(n)
	{
	}

	/*
		Copy assignment
	*/
	Matrix(const Matrix<ScalarType>& mat)
		:
		Array(mat.rows()*mat.cols(),mat.dataPtr()),
		__rows(mat.rows()),
		__cols(mat.cols())
	{
	}

	~Matrix()
	{
	}

	/*
		basic getters
	*/
	// get the number of rows
	int		rows() const {return __rows;}
	// get the number of columns
	int		cols() const {return __cols;}

	/*
		get the element of the matrix
	*/

	ScalarType& operator() (int i , int j){
		if(i<0||j<0)
			throw COException("Matrix error: index is less than zero!");
		else if (i>=__rows||j>=__cols)
			throw COException("Matrix error: index is out of range!");
		else{
			return this->__data_ptr[j*__rows+i];
		}
	}
	const ScalarType& operator() (int i,int j) const {
		return const_cast<Matrix&>(*this).operator()(i,j);
	}

	// get the element using array

	const ScalarType& data( int i ) const{
		if ( i < 0 )
			throw COException("Matrix error: index is less that zero!");
		else if ( i >= __rows*__cols )
			throw COException("Matrix error: index is out of range!");
		else
			return this->__data_ptr[i];
	}

	// set element using array

	void set ( int i , ScalarType value ){
		if ( i < 0 )
			throw COException("Matrix error: index is less that zero!");
		else if ( i >= __rows*__cols )
			throw COException("Matrix error: index is out of range!");
		else
			this->__data_ptr[i] = value;
	}

	/*
		Copy operation
	*/
	Matrix& operator= ( const Matrix<ScalarType>& mat ) {
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
		
	Matrix operator+ (const Matrix<ScalarType>& mat) {
		if ( __rows != mat.rows() || __cols != mat.cols() ) 
			throw COException("Matrix summation error: the size of two matrices are not consistent!");
		Matrix<ScalarType> result(__rows,__cols);
		for ( int i = 0 ; i < __rows*__cols ; ++ i )
			result.set(i,this->__data_ptr[i]+mat.data(i));
		return result;
	}

	// subtraction
	// need to be tested
	Matrix operator- (const Matrix<ScalarType>& mat) {
		if ( __rows != mat.rows() || __cols != mat.cols() ) 
			throw COException("Matrix subtraction error: the size of two matrices are not consistent!");
		Matrix<ScalarType> result(__rows,__cols);
		for ( int i = 0 ; i < __rows*__cols ; ++ i )
			result.set(i,this->__data_ptr[i]-mat.data(i));
		return result;
	}

	// multiply
	// need to be tested
	Vector<ScalarType> operator* ( const Vector<ScalarType>& vec ) const{
		if ( __cols != vec.size() )
			throw COException("Matrix multiply error: the size of matrix and vector are not consistent!");
		Vector<ScalarType> result(__rows);
		for ( int i = 0 ; i < __rows ; ++ i ){
			for ( int j = 0 ; j < __cols ; ++ j )
				result[i]+= operator()(i,j)*vec[j];
		}
		return result;
	}
	// need to be tested
	Matrix<ScalarType> operator* ( const Matrix<ScalarType>& mat ) const{
		if ( __cols != mat.rows() )
			throw COException("Matrix multiply error: the size of two matrices are not consistent!");
		Matrix<ScalarType> result (__rows,mat.cols());
		for ( int i = 0 ; i < __rows ; ++ i )
			for ( int j = 0 ; j < mat.cols() ; ++ j )
				for ( int k = 0 ; k < __cols ; ++ k )
					result(i,j) += operator()(i,k)*mat(k,j);
		return result;
	}
};

};
#endif