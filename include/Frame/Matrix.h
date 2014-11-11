#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.h"
namespace COPT
{
/*
	Class of 'Matrix'
		the data is stored column by column
*/
template<class FT>
class Matrix{
private:
	// the size of rows
	int			__rows;
	// the size of columns
	int			__cols;
	// the array storing the data
	FT*			__data;
public:
	// default constructor
	Matrix()
		:
		__rows(0),
		__cols(0),
		__data(NULL)
	{

	}
	Matrix(int m,int n,FT* data=NULL)
		:
		__rows(m),
		__cols(n),
		__data(new FT[m*n])
	{
		if(data){
			for ( int i = 0 ; i < m*n ; ++ i )
				__data[i] = data[i];
		}
		else{
			for ( int i = 0 ; i < m*n ; ++ i )
				__data[i] = 0.0;
		}
	}

	/*
		Copy assignment
	*/
	Matrix(const Matrix<FT>& mat)
		:
		__rows(mat.rows()),
		__cols(mat.cols()),
		__data(new FT[mat.rows()*mat.cols()])
	{
		for ( int i = 0 ; i < __rows*__cols ; ++ i )
			__data[i] = mat.data(i);
	}

	~Matrix()
	{
		SAFE_DELETE_ARRAY(__data);
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

	FT& operator() (int i , int j){
		if(i<0||j<0)
			throw COException("Matrix error: index is less than zero!");
		else if (i>=__rows||j>=__cols)
			throw COException("Matrix error: index is out of range!");
		else
			return __data[j*__cols+i];
	}
	const FT& operator() (int i,int j) const {
		return const_cast<Matrix&>(*this).operator()(i,j);
	}

	// get the element using array

	const FT& data( int i ) const{
		if ( i < 0 )
			throw COException("Matrix error: index is less that zero!");
		else if ( i >= __rows*__cols )
			throw COException("Matrix error: index is out of range!");
		else
			return __data[i];
	}

	// set element using array

	void set ( int i , FT value ){
		if ( i < 0 )
			throw COException("Matrix error: index is less that zero!");
		else if ( i >= __rows*__cols )
			throw COException("Matrix error: index is out of range!");
		else
			__data[i] = value;
	}

	/*
		Copy operation
	*/
	Matrix& operator= ( const Matrix<FT>& mat ) {
		if( __rows != mat.rows() || __cols != mat.cols() ){
			__rows = mat.rows();
			__cols = mat.cols();
			SAFE_DELETE_ARRAY(__data);
			__data = new FT[__rows*__cols];
		}

		for ( int i = 0 ; i < __rows*__cols ; ++ i )
			__data[i] = mat.data(i);
		return *this;
	}

	/*
		Mathematical operations
	*/

	// summation
	// need to be tested
		
	Matrix operator+ (const Matrix<FT>& mat) {
		if ( __rows != mat.rows() || __cols != mat.cols() ) 
			throw COException("Matrix summation error: the size of two matrices are not consistent!");
		Matrix<FT> result(__rows,__cols);
		for ( int i = 0 ; i < __rows*__cols ; ++ i )
			result.set(i,__data[i]+mat.data(i));
		return result;
	}

	// subtraction
	// need to be tested
	Matrix operator- (const Matrix<FT>& mat) {
		if ( __rows != mat.rows() || __cols != mat.cols() ) 
			throw COException("Matrix subtraction error: the size of two matrices are not consistent!");
		Matrix<FT> result(__rows,__cols);
		for ( int i = 0 ; i < __rows*__cols ; ++ i )
			result.set(i,__data[i]-mat.data(i));
		return result;
	}

	// multiply
	// need to be tested
	Vector<FT> operator* ( const Vector<FT>& vec ){
		if ( __cols != vec.size() )
			throw COException("Matrix multiply error: the size of matrix and vector are not consistent!");
		Vector<FT> result(__rows);
		for ( int i = 0 ; i < __rows ; ++ i ){
			for ( int j = 0 ; j < __cols ; ++ j )
				result[i]+= operator()(i,j)*vec[j];
		}
		return result;
	}
	// need to be tested
	Matrix<FT> operator* ( const Matrix<FT>& mat ){
		if ( __cols != mat.rows() )
			throw COException("Matrix multiply error: the size of two matrices are not consistent!");
		Matrix<FT> result (__rows,mat.cols());
		for ( int i = 0 ; i < __rows ; ++ i )
			for ( int j = 0 ; j < mat.cols() ; ++ j )
				for ( int k = 0 ; k < __cols ; ++ k )
					result(i,j) += operator()(i,k)*mat(k,j);
		return result;
	}
};

};
#endif