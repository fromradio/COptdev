// Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
// Copyright (C) MathU

#ifndef POINTER_HPP
#define POINTER_HPP

namespace COPT{
/*		Smart Pointer for COPT library
 *		A simple version of shared ptr
 */

template <class T>
class SmartPtr
{

private:
	//		the pointer to the target
	T* 			__ptr;

public:

	SmartPtr( T* ptr );


	T* operator->() const;
	T& operator*() const;

};



// Implementation of SmartPtr
template<class T>
SmartPtr<T>::SmartPtr(T * ptr)
{
	__ptr = ptr;
}

template<class T>
T* SmartPtr<T>::operator->() const 
{
	return __ptr;
}

template<class T>
T& SmartPtr<T>::operator*() const
{
	return *__ptr;
}

} // end of namespace COPT


#endif