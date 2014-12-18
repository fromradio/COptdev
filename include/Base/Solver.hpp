//		Powered by 'MathU' organization
//			Copyright@MathU

#ifndef SOLVER_H
#define SOLVER_H

/*
	A framework for a general Solver class
*/
namespace COPT{

/*		A general design for solver. The solver derives from a Time
 *		Stastistics class to help it compute time cost of the solver.
 *		The solver contains general information like max iteration 
 *		number, final iteration number, threshold, estimated error
 */	
template<class kernel , class Time=NoTimeStatistics>
class GeneralSolver
	:
	public COPTObject,
	public noncopyable
{
public:
	enum TerminalType{
		NotBeginYet,
		Optimal,
		MaxIteration,
		Unbound,
		NotFeasible
	};
protected:
	typedef typename kernel::index 				index;
	typedef typename kernel::scalar 			scalar;
	typedef typename kernel::Vector 			Vector;

	/** the journal list of solver */
	SolverJournal<GeneralSolver> 		__journal;
	/** Time computation */
	Time 								__time;
	/** max iteration number */
	index 				__max_iteration;
	/** final iteration number */
	index 				__iter_num;
	/** error threshold */
	scalar 				__thresh;
	/** final estimated error */
	scalar 				__estimated_error;
	/** the terminal type of solver */
	TerminalType 		__terminal_type;

protected:
	/** vritual function that has to be implemented by derived classes */
	/** something happens in the beginning of solving */
	virtual void solvingBegin() = 0;
	/** the real solving part */
	virtual void doSolve();
	/** the real computation part */
	virtual void doCompute() = 0;
	/** whether the terminal satisfied */
	virtual bool terminalSatisfied() const;
	
	/** an iteration 
	 *  the estimated error is returned for iteration
	 */
	virtual scalar doOneIteration() = 0;

public:
	typedef solver_object 				ObjectCategory;
	/** default constructor */
	GeneralSolver(
		const index maxiteration = 1000,
		const scalar thresh = 1e-8);
	/** virtual deconstructor */
	virtual ~GeneralSolver(){}
	/** one single iteration */
	scalar oneIteration();
	/** kernel function: solve the problem */
	void solve();
	/** kernel function: compute or initialize the problem */
	void compute();
	/** */

	/** getter and setter */
	//%{
	void 	setMaxIteration(const index maxiteration);
	index 	maxIterationNumber() const;
	void 	setIterationNumber(const index iter);
	index 	iterationNumber() const;
	void 	setThreshold( const scalar thresh);
	scalar 	threshold() const;
	void 	setEstimatedError(const scalar erro);
	scalar 	estimatedError() const;
	/** compute the objective function */
	virtual scalar objective() const = 0;
	/** get the result */
	virtual const Vector& result() const = 0;
	//%}
};

/*************Implementation of 'GeneralSolver'***************/

template<class kernel,class Time>
GeneralSolver<kernel,Time>::GeneralSolver(
	const index maxiteration,
	const scalar thresh)
	:
	__journal(*this),
	__max_iteration(maxiteration),
	__iter_num(0),
	__thresh(thresh),
	__estimated_error(0.0),
	__terminal_type(NotBeginYet)
{
}

template<class kernel,class Time>
void GeneralSolver<kernel,Time>::doSolve()
{
	__iter_num = 0;
	do
	{
		__estimated_error = this->oneIteration();
		if(++__iter_num>=__max_iteration)
		{
			__terminal_type = MaxIteration;
			break;
		}
	}while(__estimated_error>__thresh);
	if(__estimated_error<=__thresh)
		__terminal_type = Optimal;
}

template<class kernel,class Time>
bool GeneralSolver<kernel,Time>::terminalSatisfied() const
{
	return __estimated_error<__thresh;
}

template<class kernel,class Time>
void GeneralSolver<kernel,Time>::compute()
{
	__time.computationBegin();
	this->doCompute();
	__time.computationEnd();
	// after computation, the problem has not begun
	__terminal_type = NotBeginYet;
}

template<class kernel,class Time>
void GeneralSolver<kernel,Time>::solve()
{
	this->solvingBegin();
	__time.solvingBegin();
	this->doSolve();
	__time.solvingEnd();
	__journal.solveEnd();
	__time.printTimeInfo();
}

template<class kernel,class Time>
typename GeneralSolver<kernel,Time>::scalar GeneralSolver<kernel,Time>::oneIteration()
{
	__journal.iterationBegin();
	scalar e = this->doOneIteration();
	__journal.iterationEnd();
	return e;
}

template<class kernel,class Time>
void GeneralSolver<kernel,Time>::setMaxIteration( const index maxiteration)
{
	__max_iteration = maxiteration;
}

template<class kernel,class Time>
typename GeneralSolver<kernel,Time>::index GeneralSolver<kernel,Time>::maxIterationNumber() const
{
	return __max_iteration;
}

template<class kernel,class Time>
void GeneralSolver<kernel,Time>::setIterationNumber( const index iter )
{
	__iter_num = iter;
}

template<class kernel,class Time>
typename GeneralSolver<kernel,Time>::index GeneralSolver<kernel,Time>::iterationNumber() const
{
	return __iter_num;
}

template<class kernel,class Time>
void GeneralSolver<kernel,Time>::setThreshold( const scalar thresh )
{
	__thresh = thresh;
}

template<class kernel,class Time>
typename GeneralSolver<kernel,Time>::scalar GeneralSolver<kernel,Time>::threshold() const
{
	return __thresh;
}

template<class kernel,class Time>
void GeneralSolver<kernel,Time>::setEstimatedError(const scalar error )
{
	__estimated_error = error;
}

template<class kernel,class Time>
typename GeneralSolver<kernel,Time>::scalar GeneralSolver<kernel,Time>::estimatedError() const
{
	return __estimated_error;
}

//////////////////End of impelementation of 'GeneralSolver'

};

#endif