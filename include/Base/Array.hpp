// Copyright (C) Ruimin Wang ruimin.wang13@gmail.com
// Copyright (C) MathU

#ifndef ARRAY_H
#define ARRAY_H

namespace COPT
{

/*		class Array describes a base class for basic dense data types used in COPT
 *		like vector and matrix. The array can be referred to another array or independent 
 *		array. 
 */
template<class T,class S = longsize>
class Array
{
public:
	/**		define the scalar type 			*/
	typedef 				T 					ScalarType;
	/** 	define the size type 			*/
	typedef 				S 					Size;
	/**		define the category 			*/
	typedef 				data_tag 			Category;
	/** 	define the kernel trait 		*/
	typedef 				KernelTrait<T,S>	Kernel;


private:
	/** private variables */
	//%{
	/** the total size of the array */
	Size 							__size;
	/** the interval of the pointer, 1 as default */
	Size 							__inter;
	/** the pointer to the data */
	ScalarType*						__data_ptr;
	/** whether the array is referred */
	bool							__referred;
	//%}


	/* 				 protected functions
	 */
	
public:

	Array()
		:
		__size(0),
		__inter(1),
		__data_ptr(NULL),
		__referred(false)
	{
	}

	Array ( const Size size, const ScalarType* data = NULL, const Size inter = 1 )
		:
		__size(size),
		__inter(1),
		__data_ptr(new ScalarType[size]),
		__referred(false)
	{
		if ( data ){
			blas::copt_blas_copy(__size,data,1,__data_ptr,__inter);
		}
		else{
			for ( Size i = 0 ; i < __size ; ++ i )
				__data_ptr[i] = static_cast<ScalarType>(0.0);
		}
	}

	Array( const Size size , const referred_array& , ScalarType* data ,const Size inter = 1)
		:
		__size(size),
		__inter(inter),
		__data_ptr(data),
		__referred(true)
	{
	}

	Array( const Array& arr )
	{
		if(arr.isReferred())
			setReferredArray(arr.size(),arr.dataPtr(),arr.interval());
		else
			setArray(arr);
	}

	virtual ~Array()
	{
		if (__referred)
			__data_ptr = NULL;
		else
			SAFE_DELETE_ARRAY(__data_ptr);
	}


	/*
	 *		Access to data pointer for modification or other operations
	 *
	 */
	ScalarType* dataPtr(){
	 	return __data_ptr;
	}
	const ScalarType* dataPtr() const{
		return __data_ptr;
	}
	/*
	 *			Copy two arrays
	 *
	 */
	void copy( const Array& arr ) {
		resize(arr.size(),arr.interval());
		blas::copt_blas_copy(__size,arr.dataPtr(),1,__data_ptr,arr.interval());
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

	/** getter and setter */
	//%{

	/** the size of the array */ 
	const Size& size() const{return __size;}

	/** whether the array is referred */
	bool isReferred() const{return __referred;}

	/** the interval of the array */
	Size interval() const{return __inter;}

	//%}
	
	/*
	 *			resize the array to specific size
	 *
	 */
	void resize( const Size size , const Size inter = 1){
		if(__referred)
			throw COException("referred array is not allowed to be resized!");
		else{
			if(__size != size || __inter != inter ){
				__size = size;
				__inter = inter;
				SAFE_DELETE_ARRAY(__data_ptr);
				__data_ptr = new ScalarType[__size*__inter];
				for ( Size i = 0 ; i < __size ; ++ i )
					__data_ptr[i*__inter] = static_cast<ScalarType>(0.0);
			}
			else{
				for ( Size i = 0 ; i < __size ; ++ i )
					__data_ptr[i*__inter] = static_cast<ScalarType>(0.0);
			}
		}
	}

	void reset(Size size){resize(size,__inter);}

	void setArray(const Size size,const ScalarType* data,const Size inter = 1)
	{
		if ( __referred )
			throw COException("referred array is not allowed to be reset ");
		else{
			resize(size,inter);
			blas::copt_blas_copy(__size,data,1,__data_ptr,__inter);
		}
	}

	void setArray(const Array& arr)
	{
		if ( __referred )
			throw COException("referred array is not allowed to be reset ");
		else{
			resize(arr.size(),arr.interval());
			blas::copt_blas_copy(__size,arr.dataPtr(),arr.interval(),__data_ptr,__inter);
		}
	}

	/** set a referred array */
	void setReferredArray(
		const Size size,
		ScalarType* data,
		const Size inter = 1)
	{
		__referred = true;
		__size = size;
		__data_ptr = data;
		__inter = inter;
	}


	/*			Judge whether the array is valid
	 *			The array is valid if and only if the template is valid scalar type:
	 *			'float', 'double', 'std::complex<float>' or 'std::complex<double>'
	 */
	bool 			isValid() const{
		return is_scalar<ScalarType>::value;
	}

	/*
	 *
	 *
	 */
	ScalarType& operator[] ( Size i ){
		if ( i < 0 ){
			// index less than zero
			throw COException("Vector error, index less than zero.");
		}
		else if ( i >= __size ){
			// out of range
			throw COException("Vector error, index larger than the length.");
		}
		else
			return __data_ptr[i*__inter];
	}
	const ScalarType& operator[] ( Size i ) const {
		return const_cast<Array&>(*this).operator[](i);
	}

	

	/** copy assignment */
	Array& operator=(const Array& arr )
	{
		if (arr.isReferred())
		{
			__referred = true;
			__data_ptr = arr.dataPtr();
			__size = arr.size();
			__inter = arr.interval();
		}
		else{
			__referred = false;
			copy(arr);
		}
		return *this;
	}

	/**			overloaded stream */
	friend std::ostream& operator<<(std::ostream& os,const Array& arr){
		os<<"[ ";
		for ( Size i = 0 ; i< arr.size()-1 ; ++ i ){
			os<<arr[i]<<" , ";
		}
		os<<arr[arr.size()-1]<<" ]";
		return os;
	}
};


}// End of namespace COPT
#endif