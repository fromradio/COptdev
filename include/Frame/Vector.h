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
template <class T>
class Vector{
private:
	//	length of vector
	int					__length;
	//	float array
	T*					__data;
	

public:

	// the type of float number
	typedef 		T			FT;
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

	Vector( int length , FT* data = NULL )
		:
		__length(length),
		__data(new FT[length])
	{
		// __data = new FT [__length];
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
		__data(new FT[vec.size()])
	{
		for ( int i = 0 ; i < __length ; ++ i ){
			__data[i] = vec[i];
		}

	}

	/*
		API with vector in stdlib
	*/

	Vector(const std::vector<FT>& vec)
		:
		__length(static_cast<int>(vec.size())),
		__data(new FT [vec.size()])
	{
		// __data = new FT [__length];
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

	FT& operator[] ( int i ){

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
	const FT& operator[] ( int i ) const {
		return const_cast<Vector&>(*this).operator[](i);
	}
	FT& operator() ( int i ){
		return operator[](i);
	}
	const FT& operator() ( int i ) const {
		return operator[](i);
	}

	/*
		resize the Vector with specific length
			the value of each element is set as zero
	*/

	void resize( int length ){
		SAFE_DELETE_ARRAY(__data);
		__length = length;
		__data = new FT[__length];
		for ( int i = 0 ; i < __length ; ++ i) __data[i] = 0.0;
	}

	/*
		Copyt operation
	*/

	Vector& operator= (const Vector& vec ){
		__length = vec.size();
		SAFE_DELETE_ARRAY(__data);
		__data = new FT[__length];
		for ( int i = 0 ; i < __length ; ++ i ){
			__data[i] = vec[i];
		}
		return *this;
	}

	/*
		Mathematical operations
	*/

	// dot operation
	FT dot(const Vector& vec) const{
		if(__length!=vec.size()) throw COException("Vector error: the length of two vectors do not equal to each other");
		else{
			FT sum = 0;
			for ( int i = 0 ; i < __length ; ++ i ){
				sum += __data[i]*vec[i];
			}
			return sum;
		}
	}

	// summation operation
	Vector operator+ (const Vector& vec){
		if(__length!=vec.size()) throw COException("Vector error: the length of two vectors do not equal to each other");
		Vector<FT> result(__length);
		for ( int i = 0 ; i < __length ; ++ i ){
			result[i] = __data[i]+vec[i];
		}
		return result;
	}

	//subtraction operation
	Vector operator- (const Vector& vec){
		if(__length!=vec.size()) throw COException("Vector error: the length of two vectors do not equal to each other");
		Vector<FT> result(__length);
		for ( int i = 0 ; i < __length ; ++ i ){
			result[i] = __data[i]-vec[i];
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


/*
 *			compute the l2 norm of a vector
 */
// template<Vector>
// FT norm( const Vector& vec,2){
// 	FT result = 0;
// 	for ( int i = 0 ; i < vec.size() ; ++ i ){
// 		result += vec[i]*vec[i];
// 	}
// 	return sqrt(result);
// }
// template<Vector>
// FT norm( const Vector& vec,1){
// 	FT result = 0;
// 	for ( int i = 0 ; i < vec.size() ; ++ i ){
// 		result += fabs(vec[i]);
// 	}
// 	return result;
// }

};

#endif