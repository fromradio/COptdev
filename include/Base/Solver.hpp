//		Powered by 'MathU' organization
//			Copyright@MathU

#ifndef SOLVER_H
#define SOLVER_H

/*
	A framework for a general Solver class
*/
namespace COPT{


/*
	Parameter class
		the class contains basic information a solver might need
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

	/*
	 *			Copy assignment
	 *
	 */
	Parameter& operator= ( const Parameter& para){
		__initguess 	= para.initGuess();
		__maxiteration	= para.maximumIteration();
		__errorthresh	= para.errorThreshold();
	}

};


/*
 *			class for output of a solver
 */

template<class VT,class FT>
class Output
{
protected:
	// 		'true': the algorithm converges
	// 		'false': the algorithm does not converge
	bool		__converged;
	//		the result of a certain algorithm
	VT 			__result;
	//		the number of final iteration
	int 		__final_iter_num;
	//		the distance between last two iterations
	FT 			__final_error;
	//		the list that stores error in each step
	std::list<FT>
				__error_list;

public:
	Output():__converged(false),__final_iter_num(0) {}

	virtual ~Output(){}


	/*
			the getter of element of 'Output'
	*/
	//		whether the problem converges
	bool converged() {return __converged;}

	//		the result of the algorithm

	const VT&	result() {return __result;}

	//		get the list stores error

	const std::list<FT>& 	errorList() {return __error_list;}

	//		the total number of iteration

	int 	iterationNumber() {return __final_iter_num;}

	/*
	 *		basic operation:
	 *			append a new value of error
	 */
	void 	append( FT error ) {
		++ __final_iter_num;
		__error_list.push_back(error);
	}

	// algorithm ends
	//		compute the final error

	void algoEnd(){
		__final_error = __error_list.back();
	}

	// set the result

	void setResult( const VT& result ,bool converged ){
		__result = result;
		__converged = converged;
	}

	void printInfo() {
		if ( __converged )
			std::cout<<"Algorithm converges ";
		else
			std::cout<<"Algorithm fails to converge ";
	}

};

template<class VT,class FT>
class FullOutput 
	: public Output<VT,FT>
{
private:
	// a list storing every 
	std::list<VT>		__allvecs;
public:

	void append(const VT& vec,FT error){
		++ this->__final_iter_num;
		this->__error_list.push_back(error);
		__allvecs.push_back(vec);
	}
};


/*
 *	A general solver for solving problems
 *		template:
 *		VT: the input type of variable, like single variable double, vector variable Vector or matrix variable
 *		Para: the type of variable
 *		Output: the type of output containing output information
 */

// template<class kernel>
// class Solver
// 	:
// 	public COPTObject
// {
// public:
// 	typedef 	solver_tag 						Category;
// 	typedef 	typename kernel::scalar 		scalar;
// 	typedef		typename kernel::index 			index;
// 	typedef 	typename kernel::Vector 		Vector;
// 	typedef 	typename kernel::Matrix 		Matrix;

// protected:
// 	/*
// 	 *			implementation of solve
// 	 */
// 	virtual 	const Vector& doSolve (const Para& p ) = 0;
// public:
// 	// the type of variable

// 	Solver(){}

// 	virtual ~Solver(){}

// 	/*			interface of solver
// 	 *
// 	 */
// 	const Vector& solve( const Para& = Para() ){
// 		return doSolve(p);
// 	}
// };

/*
 *		an example of least square method
 *			least mean square is used
 */
// template<class VT>
// class LeastMeanSquare
// 	: public Solver<VT,Parameter<VT,typename VT::FT>, Output<VT,typename VT::FT> >
// {
// private:
// 	typedef Parameter<VT,typename VT::FT> Para;
// 	/*
// 	 *				Least square method is equal to a 
// 	 *
// 	 */
// 	virtual VT& doSolve(const Para& para){return VT();}
// public:
// };



/*
 *	an example for scalar root solver
 *	find x that satisfies that f(x)=0
 */

// template<class SFunc>
// class RootSolver:public Solver<typename SFunc::FT,Parameter<typename SFunc::FT,typename SFunc::FT>,Output<typename SFunc::FT,typename SFunc::FT> >
// {
// private:
// 	typedef				typename SFunc::FT 				FT;
// 	typedef				Parameter<FT,FT>				Para;
// 	//		reference to the target function
// 	const SFunc& 		__func;

// 	virtual const FT& doSolve(const Para& para){
// 		/*
// 			simple Netwon method
// 		*/
// 		FT xf = this->__para.initGuess();
// 		while ( fabs(__func.diff(xf))<ZERO_THRESH ){
// 			xf = xf + 0.5;
// 		}
// 		FT xn = xf-__func(xf)/__func.diff(xf);
// 		FT error = fabs(xf-xn);
// 		int iternum = 0;

// 		/*
// 		 *	the iteration terminals if the error is less than a threshold
// 		 *		or the iteration number is larger than a threshold
// 		 */
// 		while ( error > this->__para.errorThreshold() && iternum < this->__para.maximumIteration() ){
// 			xf = xn;
// 			xn = xf-__func(xf)/__func.diff(xf);
// 			error = fabs(xf-xn);
// 			++ iternum;
// 			this->__output.append(error);
// 		}
// 		this->__output.setResult(xn,true);
// 		this->__output.algoEnd();
// 		this->__output.printInfo();
// 		this->__variable = xn;

// 		return this->__variable;
// 	}
// public:
// 	//
// 	RootSolver(const SFunc& func)
// 		:
// 		__func(func)
// 	{

// 	}
// };
};

#endif