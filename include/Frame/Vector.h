//		This file is part of open library COPT
//		Copyright (c) MathU
//		Written by Ruimin Wang, ruimin.wang13@gmail.com

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
class Vector{
public:
	// the type of float number
	typedef 		FT			ScalarType;
private:
	//	length of vector
	int					__length;
	//	float array
	ScalarType*					__data;
	

public:

	// Constructor

	/*
		default constructor
	*/
	Vector()
		:
		__length(0),
		__data(NULL)
	{}

	/*
		Construct the vector with specific length
			if data is NULL, a zero vector is constructed
	*/

	Vector( int length , ScalarType* data = NULL )
		:
		__length(length),
		__data(new ScalarType[length])
	{
		// __data = new ScalarType [__length];
		if ( data ){
			for ( int i = 0 ; i < __length ; ++ i )
				__data[i] = data[i];
		}
		else {
			for ( int i = 0 ; i < __length ; ++ i )
				__data[i] = 0.0;
		}
	}

	/*
		Copy assignment
	*/
	Vector ( const Vector& vec )
		:__length(vec.size()),
		__data(new ScalarType[vec.size()])
	{
		for ( int i = 0 ; i < __length ; ++ i ){
			__data[i] = vec[i];
		}

	}

	/*
		API with vector in stdlib
	*/

	Vector(const std::vector<ScalarType>& vec)
		:
		__length(static_cast<int>(vec.size())),
		__data(new ScalarType [vec.size()])
	{
		// __data = new ScalarType [__length];
		for ( int i = 0 ; i < __length ; ++ i )
			__data[i] = vec[i];
	}

	/*
		Deconstructor
	*/

	~Vector()
	{
		SAFE_DELETE_ARRAY(__data);
	}

	/*
		get the size of the vector
	*/
		
	int size() const {
		return __length;
	}

	/*
		obtain the i-th element of the vector
	*/

	ScalarType& operator[] ( int i ){

		if ( i < 0 ){
			// index less than zero
			throw COException("Vector error, index less than zero.");
		}
		else if ( i >= __length ){
			// out of range
			throw COException("Vector error, index larger than the length.");
		}
		else
			return __data[i];
	}
	const ScalarType& operator[] ( int i ) const {
		return const_cast<Vector&>(*this).operator[](i);
	}
	ScalarType& operator() ( int i ){
		return operator[](i);
	}
	const ScalarType& operator() ( int i ) const {
		return operator[](i);
	}

	/*
		resize the Vector with specific length
			the value of each element is set as zero
	*/

	void resize( int length ){
		SAFE_DELETE_ARRAY(__data);
		__length = length;
		__data = new ScalarType[__length];
		for ( int i = 0 ; i < __length ; ++ i) __data[i] = 0.0;
	}

	/*
		Copyt operation
	*/

	Vector& operator= (const Vector& vec ){
		__length = vec.size();
		SAFE_DELETE_ARRAY(__data);
		__data = new ScalarType[__length];
		for ( int i = 0 ; i < __length ; ++ i ){
			__data[i] = vec[i];
		}
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
		for ( int i = 0 ; i < __length ; ++ i ){
			result += __data[i]*__data[i];
		}
		return result;
	}
	// dot operation
	ScalarType dot(const Vector& vec) const{
		if(__length!=vec.size()) throw COException("Vector error: the length of two vectors do not equal to each other");
		else{
			ScalarType sum = 0;
			for ( int i = 0 ; i < __length ; ++ i ){
				sum += __data[i]*vec[i];
			}
			return sum;
		}
	}

	// multiply with a scalar
	Vector operator* (ScalarType s){
		Vector result(__length);
		for ( int i = 0 ; i < __length ; ++ i ){
			result[i] = s*__data[i];
		}
		return result;
	}

	friend Vector operator* (ScalarType s,const Vector& vec){
		Vector result(vec.size());
		for ( int i = 0 ; i < vec.size() ; ++ i ){
			result[i] = s*vec[i];
		}
		return result;
	}

	// summation operation
	Vector operator+ (const Vector& vec) const{
		if(__length!=vec.size()) throw COException("Vector error: the length of two vectors do not equal to each other");
		Vector<ScalarType> result(__length);
		for ( int i = 0 ; i < __length ; ++ i ){
			result[i] = __data[i]+vec[i];
		}
		return result;
	}

	//subtraction operation
	Vector operator- (const Vector& vec) const{
		if(__length!=vec.size()) throw COException("Vector error: the length of two vectors do not equal to each other");
		Vector<ScalarType> result(__length);
		for ( int i = 0 ; i < __length ; ++ i ){
			result[i] = __data[i]-vec[i];
		}
		return result;
	}

	// 
	Vector operator- () const{
		Vector<ScalarType> result(__length);
		for ( int i = 0 ; i < __length ; ++ i ){
			result[i] = -__data[i];
		}
		return result;
	}

	// overload of stream
	friend std::ostream& operator<<(std::ostream& os,const Vector& vec){
		os<<"[ ";
		for ( int i = 0 ; i< vec.size()-1 ; ++ i ){
			os<<vec[i]<<" , ";
		}
		os<<vec[vec.size()-1]<<" ]";
		return os;
	}
};
}	// end of namespace COPT

#endif