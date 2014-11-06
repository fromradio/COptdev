// Copyright (C) Ruimin Wang ruimin.wang13@gmail.com
// Copyright (C) MathU

#ifndef ARRAY_H
#define ARRAY_H

namespace COPT
{

template<class T>
class Array
{
public:
	typedef 				T 		ScalarType;
protected:
	size_t 							__size; 		// the size of the array
	ScalarType*						__data_ptr;		// the pointer to the data itself
	
public:

	Array()
		:
		__size(0),
		__data_ptr(NULL)
	{
	}

	Array( size_t size , ScalarType* data = NULL)
		:
		__size(size),
		__data_ptr(new ScalarType[size])
	{
		if ( data ){
			for ( int i = 0 ; i < __size ; ++ i )
				__data_ptr[i] = data[i];
		}
		else{
			for ( int i = 0 ; i < __size ; ++ i )
				__data_ptr[i] = static_cast<ScalarType>(0.0);
		}
	}

	virtual ~Array()
	{
		SAFE_DELETE_ARRAY(__data_ptr);
	}


	/*
	 *		Access to data pointer for modification or other operations
	 *
	 */
	// ScalarType* dataPtr(){
	// 	return __data_ptr;
	// }
	const ScalarType* dataPtr() const{
		return __data_ptr;
	}
	/*
	 *			Copy two arrays
	 *
	 */
	void copy( const Array& arr ) {
		resize(arr.size());
		blas::copt_blas_copy(__size,arr.dataPtr(),1,__data_ptr,1);
	}

	/*
	 *			swap two arrays
	 *
	 */
	void swap ( Array& arr ) {
		if(arr.size()!=size())
			throw COException("the size of two arrays must be the same if anyone wants to swap them!");
		else
			blas::copt_blas_swap(__size,const_cast<ScalarType*>(arr.dataPtr()),1,__data_ptr,1);
	}

	/*
	 * 			the size(length) of the array
	 *
	 */
	size_t size() const{
		return __size;
	}
	/*
	 *			resize the array to specific size
	 *
	 */
	void resize( size_t size ){
		if(__size != size ){
			__size = size;
			SAFE_DELETE_ARRAY(__data_ptr);
			__data_ptr = new ScalarType[__size];
			for ( int i = 0 ; i < __size ; ++ i )
				__data_ptr[i] = static_cast<ScalarType>(0.0);
		}
	}
	/*			Judge whether the array is valid
	 *			The array is valid if and only if the template is valid scalar type:
	 *			'float', 'double', 'std::complex<float>' or 'std::complex<double>'
	 */
	bool 			isValid(){
		return {is_scalar<ScalarType>::value;}
	}
};


}// End of namespace COPT
#endif