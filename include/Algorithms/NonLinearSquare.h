#ifndef NON_LINEAR_SOLVER_H
#define NON_LINEAR_SOLVER_H

#include "Solver.h"


namespace COPT
{
typedef Vector<double>													VT;
typedef Solver<VT,Parameter<VT,double>,Output<VT,double> >				Sol;

class NonLinearSquare
	: public Sol
{
private:
	/*
	 *				the list of vector functions that are used
	 *				
	 */
	std::list<VectorFunction<VT> *>				__vfs;

	/*
	 *				non-linear square solver
	 */
	virtual VT& doSolve(const Para& para)
	{
		return VT();
	}
public:
	/*
	 *				Constructor
	 */
	NonlinearSquare();
};
};

#endif