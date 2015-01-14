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

#ifndef ARRAY_IMPL_HPP__
#define ARRAY_IMPL_HPP__

namespace COPT
{

/*********************Implementation of 'Array'************************/

template<class scalar,class index>
Array<scalar,index>::Array()
	:
	__size(0),
	__inter(1),
	__data_ptr(nullptr),
	__referred(false)
{
}

template<class scalar,class index>
Array<scalar,index>::Array(const index size, const scalar* data, const index inter)
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

template<class scalar,class index>
Array<scalar,index>::Array(const index size, const referred_array&, scalar *data, const index inter)
	:
	__size(size),
	__inter(inter),
	__data_ptr(data),
	__referred(true)
{
}

template<class scalar,class index>
Array<scalar,index>::Array(const Array& arr)
{
	if(arr.isReferred())
		setReferredArray(arr.size(),arr.dataPtr(),arr.interval());
	else
		setArray(arr);
}

template<class scalar,class index>
Array<scalar,index>::~Array()
{
	if (__referred)
		__data_ptr = NULL;
	else
		SAFE_DELETE_ARRAY(__data_ptr);
}

template<class scalar,class index>
void Array<scalar,index>::clear()
{
	if (__referred)
	{
		__data_ptr = NULL;
		__size = 0;
		__referred = false;
	}
	else
	{
		SAFE_DELETE_ARRAY(__data_ptr);
		__size = 0;
	}
}

template<class scalar,class index>
scalar* Array<scalar,index>::dataPtr()
{
	return __data_ptr;
}

template<class scalar,class index>
const scalar* Array<scalar,index>::dataPtr() const
{
	return __data_ptr;
}

template<class scalar,class index>
void Array<scalar,index>::copy(const Array& arr)
{
	if(this->isReferred())
		__referred = false;
	resize(arr.size(),arr.interval());
	blas::copt_blas_copy(__size,arr.dataPtr(),1,__data_ptr,arr.interval());
}

template<class scalar,class index>
void Array<scalar,index>::swap(Array& arr)
{
	if(arr.size()!=size())
		throw COException("the size of two arrays must be the same if anyone wants to swap them!");
	else if(arr.isReferred()||this->isReferred())
		throw COException("Swap requires that both arrays are not referred array");
	else
		blas::copt_blas_swap(__size,const_cast<scalar*>(arr.dataPtr()),1,__data_ptr,1);
}

template<class scalar,class index>
const index& Array<scalar,index>::size() const
{
	return __size;
}

template<class scalar,class index>
bool Array<scalar,index>::isReferred() const
{
	return __referred;
}

template<class scalar,class index>
index Array<scalar,index>::interval() const
{
	return __inter;
}

template<class scalar,class index>
void Array<scalar,index>::resize(const index size, const index inter)
{
	if(__referred)
		throw COException("referred array is not allowed to be resized!");
	else
	{
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

template<class scalar,class index>
void Array<scalar,index>::reset(const index size,const index inter)
{
	if(this->isReferred())
	{
		clear();
		__data_ptr = new scalar[__size*__inter];
		for ( index i = 0 ; i < __size ; ++ i )
			__data_ptr[i*__inter] = static_cast<scalar>(0.0);
	}
	else
	{
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

template<class scalar,class index>
void Array<scalar,index>::setArray(const index size, const scalar *data, const index inter)
{
	if (__referred)
		throw COException("referred array is not allowed to be reset ");
	else{
		resize(size,inter);
		blas::copt_blas_copy(__size,data,1,__data_ptr,__inter);
	}
}

template<class scalar,class index>
void Array<scalar,index>::setArray(const Array &arr)
{
	if ( __referred )
		throw COException("referred array is not allowed to be reset ");
	else{
		resize(arr.size(),arr.interval());
		blas::copt_blas_copy(__size,arr.dataPtr(),arr.interval(),__data_ptr,__inter);
	}
}

template<class scalar,class index>
void Array<scalar,index>::setReferredArray(const index size, scalar* data, const index inter)
{
	__referred = true;
	__size = size;
	__data_ptr = data;
	__inter = inter;
}

template<class scalar,class index>
bool Array<scalar,index>::isValid() const
{
	return is_scalar<scalar>::value;
}

template<class scalar,class index>
scalar& Array<scalar,index>::operator[](index i)
{
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

template<class scalar,class index>
const scalar& Array<scalar,index>::operator[](index i)const
{
	return const_cast<Array&>(*this).operator[](i);
}

template<class scalar,class index>
Array<scalar,index>& Array<scalar,index>::operator=(const Array& arr)
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

}

#endif