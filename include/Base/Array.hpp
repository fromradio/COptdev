// This file is part of COPT, a light-weight C++ based optimization open source library
//
// Copyright (C) 2015 Ruimin Wang <ruimin.wang13@gmail.com>
// Copyright (C) 2015 MathU
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#ifndef ARRAY_HPP__
#define ARRAY_HPP__

namespace COPT
{
/*		
 *
 *
 */
class DataObject: public COPTObject
{

};
/*		class Array describes a base class for basic dense data types used in COPT
 *		like vector and matrix. The array can be referred to another array or independent 
 *		array. 
 */

template<class T,class I>
class Array
	:
	public DataObject
{
public:
	/**		define the scalar type 			*/
	typedef 				T 					scalar;
	/** 	define the size type 			*/
	typedef 				I 					index;
	/**		define the category 			*/
	typedef 				array_object 		ObjectCategory;
	/** 	define the kernel trait 		*/
	typedef 				KernelTrait<T,I>	Kernel;


private:
	/** private variables */
	//%{
	/** the total size of the array */
	index 							__size;
	/** the interval of the pointer, 1 as default */
	index 							__inter;
	/** the pointer to the data */
	scalar*							__data_ptr;
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

	Array ( const index size, const scalar* data = NULL, const index inter = 1 )
		:
		__size(size),
		__inter(1),
		__data_ptr(new scalar[size]),
		__referred(false)
	{
		if ( data ){
			blas::copt_blas_copy(__size,data,1,__data_ptr,__inter);
		}
		else{
			for ( index i = 0 ; i < __size ; ++ i )
				__data_ptr[i] = static_cast<scalar>(0.0);
		}
	}

	Array( const index size , const referred_array& , scalar* data ,const index inter = 1)
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
	scalar* dataPtr(){
	 	return __data_ptr;
	}
	const scalar* dataPtr() const{
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
			blas::copt_blas_swap(__size,const_cast<scalar*>(arr.dataPtr()),1,__data_ptr,1);
	}

	/** getter and setter */
	//%{

	/** the size of the array */ 
	const index& size() const{return __size;}

	/** whether the array is referred */
	bool isReferred() const{return __referred;}

	/** the interval of the array */
	index interval() const{return __inter;}

	//%}
	
	/*
	 *			resize the array to specific size
	 *
	 */
	void resize( const index size , const index inter = 1){
		if(__referred)
			throw COException("referred array is not allowed to be resized!");
		else{
			if(__size != size || __inter != inter ){
				__size = size;
				__inter = inter;
				SAFE_DELETE_ARRAY(__data_ptr);
				__data_ptr = new scalar[__size*__inter];
				for ( index i = 0 ; i < __size ; ++ i )
					__data_ptr[i*__inter] = static_cast<scalar>(0.0);
			}
			else{
				for ( index i = 0 ; i < __size ; ++ i )
					__data_ptr[i*__inter] = static_cast<scalar>(0.0);
			}
		}
	}

	void reset(index size){resize(size,__inter);}

	void setArray(const index size,const scalar* data,const index inter = 1)
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
		const index size,
		scalar* data,
		const index inter = 1)
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
		return is_scalar<scalar>::value;
	}

	/*
	 *
	 *
	 */
	scalar& operator[] ( index i ){
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
	const scalar& operator[] ( index i ) const {
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
		for ( index i = 0 ; i< arr.size()-1 ; ++ i ){
			os<<arr[i]<<" , ";
		}
		os<<arr[arr.size()-1]<<" ]";
		return os;
	}
};


}// End of namespace COPT
#endif