//		This file is part of open library COPT
//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <iostream>

/*
	basic operation that will be used in the pipeline
*/
#include "BasicOperation.h"

/*
	the exception class for COPT
*/
#include "COException.h"

namespace COPT
{
template <class FT>
class Vector : public Array<FT>{
public:
	// the type of float number
	typedef 				Array<FT>				Array;
	typedef typename 		Array::ScalarType		ScalarType;

public:

	// Constructor

	/*
		default constructor
	*/
	Vector()
		:
		Array()
	{}

	/*
		Construct the vector with specific length
			if data is NULL, a zero vector is constructed
	*/

	Vector( int length , ScalarType* data = NULL )
		:
		Array(length,data)
	{
	}

	/*
		Copy assignment
	*/
	Vector ( const Vector& vec )
		
	{
		this->__size=vec.size();
		this->__data_ptr=new ScalarType[vec.size()];
		blas::copt_blas_copy(this->__size,vec.dataPtr(),1,this->__data_ptr,1);
	}

	/*
		API with vector in stdlib
	*/

	Vector(const std::vector<ScalarType>& vec)
	{
		this->__size=static_cast<size_t>(vec.size());
		this->__data_ptr=new ScalarType [vec.size()];
		for ( int i = 0 ; i < this->__size ; ++ i )
			this->__data_ptr[i] = vec[i];
	}

#ifdef EIGEN
	Vector(const Eigen::Matrix<ScalarType,Eigen::Dynamic,1>& vec)
		:
		Array(vec.size())
	{
		for ( int i = 0 ; i < this->__size ; ++ i ){
			this->__data_ptr[i]  =vec(i);
		}
	}
#endif

	/*
		Deconstructor
	*/

	~Vector()
	{
	}

	/*
		Copyt operation
	*/

	Vector& operator= (const Vector& vec ){
		this->resize(vec.size());
		this->setData(vec.size(),vec.dataPtr());
		return *this;
	}

	/*
		Mathematical operations
	*/
	/*
	 * 			Square norm of the vector
	 */
	ScalarType squaredNorm() const{
		ScalarType result = 0;
		for ( int i = 0 ; i < this->__size ; ++ i ){
			result += this->__data_ptr[i]*this->__data_ptr[i];
		}
		return result;
	}
	// dot operation
	ScalarType dot(const Vector& vec) const{
		if(this->__size!=vec.size()) throw COException("Vector error: the length of two vectors do not equal to each other");
		else{
			ScalarType sum = blas::copt_blas_dot(this->__size,this->__data_ptr,1,vec.dataPtr(),1);
			return sum;
		}
	}

	// scale with a special length
	void scale(ScalarType s){
		blas::copt_blas_scal(this->__size,s,this->__data_ptr,1);
	}
	// multiply with a scalar
	Vector operator* (ScalarType s){
		Vector result(*this);
		result.scale(s);
		return result;
	}

	friend Vector operator* (ScalarType s,const Vector& vec){
		Vector result(vec);
		result.scale(s);
		return result;
	}

	// summation operation
	Vector operator+ (const Vector& vec) const{
		if(this->__size!=vec.size()) throw COException("Vector error: the length of two vectors do not equal to each other");
		Vector<ScalarType> result(this->__size);
		for ( int i = 0 ; i < this->__size ; ++ i ){
			result[i] = this->__data_ptr[i]+vec[i];
		}
		return result;
	}

	//subtraction operation
	Vector operator- (const Vector& vec) const{
		if(this->__size!=vec.size()) throw COException("Vector error: the length of two vectors do not equal to each other");
		Vector<ScalarType> result(this->__size);
		for ( int i = 0 ; i < this->__size ; ++ i ){
			result[i] = this->__data_ptr[i]-vec[i];
		}
		return result;
	}

	// 
	Vector operator- () const{
		Vector<ScalarType> result(this->__size);
		for ( int i = 0 ; i < this->__size ; ++ i ){
			result[i] = -this->__data_ptr[i];
		}
		return result;
	}

	/* overload of stream
	 */
	friend std::ostream& operator<<(std::ostream& os,const Vector& vec){
		os<<"[ ";
		for ( int i = 0 ; i< vec.size()-1 ; ++ i ){
			os<<vec[i]<<" , ";
		}
		os<<vec[vec.size()-1]<<" ]";
		return os;
	}


	/*		Generate special vectors
	 */
	/*			Generate e vector = [0,0,...,1,0,...]
	 *			/param size:		the size of the vector
	 *			/param i:			the index of non-zero element
	 */
	static Vector vecE(size_t size,int i)
	{
		if ( i < 0 || i >= size ) throw COException("Index error: out of range!");
		Vector vec(size);
		vec[i] = 1.0;
		return vec;
	}
	/*			Generate e vector = [0,0,...,s,0,...]
	 *			/param size:		the size of the vector
	 *			/param i:			the index of non-zero element
	 *			/param s:			the value of non-zero element
	 */
	static Vector vecE(size_t size,int i,const ScalarType s)
	{
		if ( i < 0 || i >= size ) throw COException("Index error: out of range!");
		Vector vec(size);
		vec[i] = s;
		return vec;
	}
};

}	// end of namespace COPT

#endif