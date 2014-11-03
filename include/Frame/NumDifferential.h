#ifndef NUM_DIFFERENTIAL_H
#define NUM_DIFFERENTIAL_H

#include "Matrix.h"
/*
	compute the differential of functions
*/
namespace COPT{
template<class SFunc>
class ScalarDifferential{
private:
	// the type of float number
	typedef typename SFunc::FT 			FT;
	// the reference to the scalar function
	const SFunc&						__func;
	// the error threshold
	FT 									__epsilon;


	// the iteration number that is used
	// this parameter is only used in dev version
	int 								__iterused;

	/*
		private Functions
	*/
	// central differentce
	FT centerDifference (FT x,FT h){
		return (__func(x+h)-__func(x-h))/(2*h);
	}
public:
	/*
		Constructor:
			no default constructor is allowed
			one must specify the function that is used
	*/
	ScalarDifferential(const SFunc& func,FT epsilon = 1e-5)
		:
		__func(func),
		__epsilon(epsilon),
		__iterused(0)
	{}

	// compute the differential
	FT diff(FT x, FT h = 0.01){
		FT error 		= __epsilon+1.0;
		FT forwarddiff 	= centerDifference(x,h);
		FT currdiff		= forwarddiff;
		__iterused		= 0;
		while ( error > __epsilon ){
			h /= 2.0;
			forwarddiff = currdiff;
			currdiff 	= centerDifference(x,h);
			error 		= fabs(currdiff - forwarddiff);
			if ( ++ __iterused > 20 ){
				std::cerr<<"Reach the max iteration number, differential result might be inaccurate"<<std::endl;
				break;
			}
		}
		return currdiff;
	}

};
};

#endif