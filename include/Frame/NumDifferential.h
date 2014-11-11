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

	// second order differential
	FT secondDifference (FT x,FT h){
		return (__func(x+h)+__func(x-h)-2*__func(x))/(h*h);
	}

public:
	/*
	 *	Constructor:
	 *		no default constructor is allowed
	 *		one must specify the function that is used
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

	// compute second order differential
	FT sDiff(FT x,FT h = 0.01){
		FT error 			= __epsilon + 1.0;
		FT forwarddiff 		= secondDifference(x,h);
		FT currdiff 		= forwarddiff;
		__iterused 			= 0;
		while (error > __epsilon ){
			h /= 2.0;
			forwarddiff = currdiff;
			currdiff 	= secondDifference(x,h);
			error 		= fabs(currdiff - forwarddiff);
			if ( ++ __iterused > 20 ){
				std::cerr<<"Reach the max iteration number, differential result might be inaccurate"<<std::endl;
				break;
			}
		}
		return currdiff;
	}

};


template<class VFunc>
class VectorDifferential{
private:
	typedef typename VFunc::Vector 					Vector;
	typedef typename Vector::FT 					FT;
	/*
	 *		const reference to the function
	 */
	const VFunc& 									__vfunc;
	// 
	FT 												__epsilon;
	//
	int 											__iterused;

	/*
	 *				Private Function
	 */
	FT centerDifference(const Vector& vec,FT h,int i){
		Vector vp(vec),vm(vec);
		vp[i] += h;
		vm[i] -= h;
		return (__vfunc(vp)-__vfunc(vm))/(2*h);
	}

	/*
	 *
	 *
	 *
	 */
	FT computeDiff(const Vector& vec,int i,FT h = 0.01){
		FT error 		= __epsilon+1.0;
		FT forwarddiff 	= centerDifference(vec,h,i);
		FT currdiff		= forwarddiff;
		__iterused		= 0;
		while ( error > __epsilon ){
			h /= 2.0;
			forwarddiff = currdiff;
			currdiff 	= centerDifference(vec,h,i);
			error 		= fabs(currdiff - forwarddiff);
			if ( ++ __iterused > 20 ){
				std::cerr<<"Reach the max iteration number, differential result might be inaccurate"<<std::endl;
				break;
			}
		}
		return currdiff;
	}

	/*
	 *				compute differential of second order
	 */
	// FT secondDifference(const Vector& vec,FT h,int i,int j){
	// 	Vector vp(vec),vm(vec);
	// 	vp[i] +
	// }
public:
	VectorDifferential(const VFunc& func,FT epsilon = 1e-5) 
		: 
		__vfunc(func),
		__epsilon(epsilon)
		{}
	/*
	 *			Main function, compute the gradient of '__func'
	 */
	Vector gradient(const Vector& vec,FT h = 0.01){
	 	Vector result(vec.size());
	 	for ( int i = 0 ; i < vec.size() ; ++ i ){
	 		result[i] = computeDiff(vec,i,h);
	 	}
	 	return result;
	 }
};
};

#endif