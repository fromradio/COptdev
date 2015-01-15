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
	
public:

	/** constructors and deconstructor */
	//%{
	/** default constructor */
	Array();
	/** constructor for a standard array. */
	Array (const index size, const scalar *data = nullptr, const index inter = 1);
	/** constructor for a referred array. */
	Array(const index size, const referred_array&, scalar* data, const index inter = 1);
	/** copy constructor */
	Array(const Array& arr);
	/** deconstructor */
	virtual ~Array();
	//%}

	/** clear the array */
	void clear();

	/**	Access to data pointer for modification or other operations */
	scalar* dataPtr();
	const scalar* dataPtr() const;

	/**	copy from another array, 
	  * unlike copy constructor, the array will not be a reffered one
	   */
	void copy(const Array& arr);

	/**	swap two arrays */
	void swap (Array& arr);


	/** getter and setter */
	//%{
	/** the size of the array */ 
	const index& size() const;
	/** whether the array is referred */
	bool isReferred() const;
	/** the interval of the array */
	index interval() const;
	//%}
	
	/**		resize the array to specific size
	 *
	 */
	void resize(const index size, const index inter = 1);

	/** reset the array even if the array is a referred array */
	void reset(const index size, const index inter=1);

	/** set the array with given size and data */
	void setArray(const index size, const scalar* data, const index inter = 1);

	/** set the array with given array */
	void setArray(const Array& arr);

	/** set a referred array */
	void setReferredArray(
		const index size,
		scalar* data,
		const index inter = 1);


	/*	Judge whether the array is valid
	 *	The array is valid if and only if the template is valid scalar type:
	 *	'float', 'double', 'std::complex<float>' or 'std::complex<double>'
	 */
	bool isValid() const;

	/** data access*/
	scalar& operator[] (index i);
	const scalar& operator[] (index i)const;

	/** copy assignment */
	Array& operator=(const Array& arr);

	/** overloaded stream */
	friend std::ostream& operator<<(std::ostream& os, const Array& arr){
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