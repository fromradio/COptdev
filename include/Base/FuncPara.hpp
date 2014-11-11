#ifndef FUNC_PARA_H
#define FUNC_PARA_H

#include <iostream>

template<class Scalar>
class FuncPara{
public:
	//
	FuncPara 		()
		:
		__length(0),
		__data(NULL)
	{
	}


	FuncPara 		(int length)
		:__length(length),
		__data(NULL)
	{
		resize(length);
	}

	~FuncPara		()
	{
		if(__data) delete __data;
	}

	/*
		setter and getter
	*/
	int  size	() 	const
	{
		return __length;
	}


	// get the i-th number
	Scalar& 		operator[]	( int i )
	{
		if(i>=0&&i<__length)
			return __data[i];
		else if(i<0){
			std::cerr<<"Error: The index is less than zero, the first value is returned"<<std::endl;
			return __data[0];
		}
		else{
			std::cerr<<"Error: the index is out of range, the last value is returned"<<std::endl;
			return __data[__length-1];
		}	
	}


	const Scalar& 	operator[]	( int i )	const
	{
		if(i>=0&&i<__length)
			return __data[i];
		else if(i<0){
			std::cerr<<"Error: The index is less than zero, the first value is returned"<<std::endl;
			return __data[0];
		}
		else{
			std::cerr<<"Error: the index is out of range, the last value is returned"<<std::endl;
			return __data[__length-1];
		}
	}


	// resize the vector
	void 			resize		( int length )
	{
		if ( length < 0 ) {
			std::cerr<<"Error: the size of para is less than zero"<<std::endl;
			return;
		}
		__length = length;
		if (__data) delete __data;
		__data = new Scalar(__length);
	}

	// set zero to all elements
	void 			setZero		(  )
	{
		for ( int i = 0 ; i < __length ; ++ i )
		{
			__data[i] = 0.0;
		}
	}


	// set according to a vector
	void	setVector	( int length , Scalar* vec)
	{
		if (__data) delete __data;
		__length = length;
		__data = new Scalar(__length);
		for ( int i = 0 ; i < __length ; ++ i )
			__data[i] = vec[i];
	}
private:
	int 			__length;
	Scalar*			__data;
};



#endif