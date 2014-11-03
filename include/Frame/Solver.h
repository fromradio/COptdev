//		Powered by 'MathU' organization
//			Copyright@MathU

#ifndef SOLVER_H
#define SOLVER_H

#include "BasicOperation.h"
/*
	A framework for a general Solver class
*/
namespace COPT{


/*
	Parameter class
		the class contains basic information a solver needs
*/

template<class FT,class VT = FT>
class Parameter{
private:
	//		the initial guess of a problem
	VT 				__initguess;
	//		maximum iteration number
	int 			__maxiteration;
	//		error threshold
	FT				__errorthresh;
	//		the dimension of the problem
	int 			__dim;
	//		whether initial guess is given
	bool			__giveninit;
public:
	//		Default constructor
	Parameter()
		: 
		__maxiteration(100),
		__errorthresh(1e-5),
		__giveninit(false)
	{}

	//		Constructor in which initial guess is given

	Parameter( const VT& initial, int maxiteration = 100, FT errorthresh = 1e-5 )
		:
		__initguess(initial),
		__maxiteration(maxiteration),
		__errorthresh(errorthresh),
		__giveninit(true)
	{}

	//		Deconstructor

	~Parameter(){}

	// get the init guess
	const VT& 		initGuess()const {
		if ( __giveninit )
			return __initguess;
		else{
			if (typeid(FT)==typeid(VT))
				return 0.0;
			else
				return VT();
		}
	}

	// get the maximum iteration number
	const int&		maximumIteration() const {
		return __maxiteration;
	}

	// get the error threshold
	const FT&		errorThreshold() const {
		return __errorthresh;
	}

};


/*
	class for output of a solver
*/

template<class VT,class FT>
class Output
{
private:
	// 		'true': the algorithm converges
	// 		'false': the algorithm does not converge
	bool		__converged;
	//		the result of a certain algorithm
	VT 			__result;
public:
};



/*
	A general solver for solving problems
		template:
		VT: the input type of variable, like single variable double, vector variable Vector or matrix variable
		Para: the type of variable
		Output: the type of output containing output information
*/

template<class VT,class Para,class Output>
class Solver{
protected:
	//	variable
	VT 			__variable;
	// parameter
	Para		__para;
	// output
	Output 		__output;
public:
	// the type of variable
	typedef			VT 			VaType;

	Solver(){}
	virtual ~Solver(){}
	virtual const VT& solve( const Para& p = Para()) = 0;
};


/*
	an example for scalar root solver
	find x that satisfies that f(x)=0
*/

template<class SFunc>
class RootSolver:public Solver<typename SFunc::FT,Parameter<typename SFunc::FT,typename SFunc::FT>,Output<typename SFunc::FT,typename SFunc::FT> >
{
private:
	typedef				typename SFunc::FT 				FT;
	typedef				Parameter<FT,FT>				Para;
	//		reference to the target function
	const SFunc& 		__func;
public:
	//
	RootSolver(const SFunc& func)
		:
		__func(func)
	{
	}

	//			compute the root of SFunc
	const FT& solve( const Para& p = Para() ){
		/*
			Netwon method
		*/
		FT xf = this->__para.initGuess();
		while ( fabs(__func.diff(xf))<ZERO_THRESH ){
			xf = xf + 0.5;
		}
		FT xn = xf-__func(xf)/__func.diff(xf);
		FT error = fabs(xf-xn);
		int iternum = 0;
		while ( error > this->__para.errorThreshold() && iternum < this->__para.maximumIteration() ){
			xf = xn;
			xn = xf-__func(xf)/__func.diff(xf);
			error = fabs(xf-xn);
			++ iternum;
		}
		this->__variable = xn;
		return this->__variable;
	}
};
};

#endif