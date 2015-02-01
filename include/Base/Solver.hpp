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


#ifndef SOLVER_HPP__
#define SOLVER_HPP__

/*
	A framework for a general Solver class
*/
namespace COPT{

/** 		Description for objective function. 
  * 		Simple use case:
  *			double func(const Matrix& mat){return mat.frobeniusNorm();}
  *			ObjectiveFunction f = func;
  *			Matrix m = Matrix::identity(4,4);
  *			std::cout<<f(m); // 2.0;
  */
template<class Scalar,class ArgType>
using ObjectiveFunction = Scalar(*)(const ArgType&);

template<class Scalar,class ArgType,class Option,class Para>
using OneIteration = Scalar(*)(ArgType&,Option&,Para&);

template<class ArgType,class Option>
using CheckTermination = bool(*)(const ArgType&,const Option&);

// template<class Scalar,class ArgType,ObjectiveFunction<Scalar,ArgType> ObFunc>
// class Solver
// {
// private:
// 	/** the result */
// 	ArgType 		__x;
// public:
// 	Scalar objective() const;
// 	static Scalar objective(const ArgType& x);
// };

// template<class Scalar,class ArgType,ObjectiveFunction<Scalar,ArgType> ObFunc>
// Scalar Solver<Scalar,ArgType,ObFunc>::objective()const
// {
// 	return ObFunc(__x);
// }

// template<class Scalar,class ArgType,ObjectiveFunction<Scalar,ArgType> ObFunc>
// Scalar Solver<Scalar,ArgType,ObFunc>::objective(const ArgType& x)
// {
// 	return ObFunc(x);
// }

/** types of final termination */
enum TerminalType
{
	C, 		// converged
	U,		// Unbound
	N,		// not feasible
	M 		// maximum iteration is reached
};

template<class Scalar>
struct BasicOption
{
	/** the maximum iteration number */
	int 					MaxIter;
	/** the error threshold */
	Scalar 					Threshold;
};

template<class Scalar>
struct BasicParameter
{
	/** the iteration number */
	int 					IterNum;
	/** the iterative error */
	Scalar 					Error;
	/** the computation time */
	double 					ComputationTime;
	/** the terminal type */
	TerminalType 			Termination;
};



template<class Scalar,class ArgType,class OutputType=ArgType,class Option=BasicOption<Scalar>,class Parameter = BasicParameter<Scalar> >
class Solver
{
private:

	typedef ObjectiveFunction<Scalar,ArgType> 		ObjectiveFunction;
	typedef OneIteration<Scalar,ArgType,Option,Parameter> 			IterationFunction;
	typedef CheckTermination
	

	ArgType 				__x;
	/** the result */
	OutputType 				__result;
	/** the option of the solver */
	Option 					__op;
	/** the parameter */
	Parameter 				__para;
	/** the call for objective function */
	ObjectiveFunction 		__ob_func;
	/** the iteration for the solver */
	IterationFunction 		__iter_func;

	virtual void doSolve();

public:

	Solver(ObjectiveFunction obfunc = nullptr,IterationFunction iterfunc = nullptr);

	/** return the result of the solver */
	OutputType result() const;


	/** set the objective function */
	void setObjectiveFunction(ObjectiveFunction func);
	/** return the current objective value */
	Scalar objective() const;
	/** how to compute objective function */
	Scalar objective(const ArgType& x) const;

	/** set the iteration function */
	void setIterationFunction(IterationFunction func);

	/** get the option */
	void Option& option() const;

	/** get the parameter */
	const Parameter& parameter() const;
};

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
Solver<Scalar,ArgType,OutputType,Option,Parameter>::Solver(ObjectiveFunction obfunc,IterationFunction iterfunc)
	:
	__ob_func(obfunc),
	__iter_func(iterfunc)
{
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
OutputType Solver<Scalar,ArgType,OutputType,Option,Parameter>::result() const
{
	return __result;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
void Solver<Scalar,ArgType,OutputType,Option,Parameter>::setObjectiveFunction(ObjectiveFunction func)
{
	__ob_func = func;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
Scalar Solver<Scalar,ArgType,OutputType,Option,Parameter>::objective() const
{
	if (__ob_func)
		return __ob_func(__x);
	else
		throw COException("Solver error: Objective Function is not defined yet!");
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
Scalar Solver<Scalar,ArgType,OutputType,Option,Parameter>::objective(const ArgType& x) const
{
	if (__ob_func)
		return __ob_func(x);
	else
		throw COException("Solver error: Objective Function is not defined yet!");
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
void Solver<Scalar,ArgType,OutputType,Option,Parameter>::setIterationFunction(IterationFunction func)
{
	__iter_func = func;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
void Solver<Scalar,ArgType,OutputType,Option,Parameter>::doSolve()
{
	clock_t b,e;
	b = clock();
	for(int i=0; i<__op.MaxIter; ++i)
	{
		__iter_func(__x,__op,__para);
	}
	e = clock();
	if(i==__op.MaxIter) 
		__para.Termination = M; // max iteration number

}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
const Option& Solver<Scalar,ArgType,OutputType,Option,Parameter>::option() const
{
	return __op;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
const Parameter& Solver<Scalar,ArgType,OutputType,Option,Parameter>::parameter() const
{
	return __para;
}


/*			A general design for solver. The solver derives from a Time
 *			Stastistics class to help it compute time cost of the solver.
 *			The solver contains general information like max iteration 
 *			number, final iteration number, threshold, estimated error
 */	
template<class kernel,class Time=NoTimeStatistics,class ResultType=typename kernel::Vector>
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
	typedef typename kernel::podscalar 			podscalar;
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
	podscalar			__thresh;
	/** final estimated error */
	podscalar 			__estimated_error;
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
	virtual podscalar doOneIteration() = 0;

public:
	typedef solver_object 				ObjectCategory;
	/** default constructor */
	GeneralSolver(
		const index maxiteration = 1000,
		const podscalar thresh = 1e-8);
	/** virtual deconstructor */
	virtual ~GeneralSolver(){}
	/** one single iteration */
	podscalar oneIteration();
	/** kernel function: solve the problem */
	void solve();
	/** kernel function: compute or initialize the problem */
	void compute();
	/** */

	/** getter and setter */
	//%{
	void 		setMaxIteration(const index maxiteration);
	index 		maxIterationNumber() const;
	void 		setIterationNumber(const index iter);
	index 		iterationNumber() const;
	void 		setThreshold( const podscalar thresh);
	podscalar 	threshold() const;
	void 		setEstimatedError(const podscalar error);
	podscalar 	estimatedError() const;
	/** compute the objective function */
	virtual podscalar objective() const = 0;
	/** get the result */
	virtual const ResultType& result() const = 0;
	//%}
};

/*************Implementation of 'GeneralSolver'***************/

template<class kernel,class Time,class ResultType>
GeneralSolver<kernel,Time,ResultType>::GeneralSolver(
	const index maxiteration,
	const podscalar thresh)
	:
	__journal(*this),
	__max_iteration(maxiteration),
	__iter_num(0),
	__thresh(thresh),
	__estimated_error(0.0),
	__terminal_type(NotBeginYet)
{
}

template<class kernel,class Time,class ResultType>
void GeneralSolver<kernel,Time,ResultType>::doSolve()
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

template<class kernel,class Time,class ResultType>
bool GeneralSolver<kernel,Time,ResultType>::terminalSatisfied() const
{
	return __estimated_error<__thresh;
}

template<class kernel,class Time,class ResultType>
void GeneralSolver<kernel,Time,ResultType>::compute()
{
	__time.computationBegin();
	this->doCompute();
	__time.computationEnd();
	// after computation, the problem has not begun
	__terminal_type = NotBeginYet;
}

template<class kernel,class Time,class ResultType>
void GeneralSolver<kernel,Time,ResultType>::solve()
{
	this->solvingBegin();
	__time.solvingBegin();
	this->doSolve();
	__time.solvingEnd();
	__journal.solveEnd();
	__time.printTimeInfo();
}

template<class kernel,class Time,class ResultType>
typename GeneralSolver<kernel,Time,ResultType>::podscalar GeneralSolver<kernel,Time,ResultType>::oneIteration()
{
	__journal.iterationBegin();
	podscalar e = this->doOneIteration();
	__journal.iterationEnd();
	return e;
}

template<class kernel,class Time,class ResultType>
void GeneralSolver<kernel,Time,ResultType>::setMaxIteration( const index maxiteration)
{
	__max_iteration = maxiteration;
}

template<class kernel,class Time,class ResultType>
typename GeneralSolver<kernel,Time,ResultType>::index GeneralSolver<kernel,Time,ResultType>::maxIterationNumber() const
{
	return __max_iteration;
}

template<class kernel,class Time,class ResultType>
void GeneralSolver<kernel,Time,ResultType>::setIterationNumber( const index iter )
{
	__iter_num = iter;
}

template<class kernel,class Time,class ResultType>
typename GeneralSolver<kernel,Time,ResultType>::index GeneralSolver<kernel,Time,ResultType>::iterationNumber() const
{
	return __iter_num;
}

template<class kernel,class Time,class ResultType>
void GeneralSolver<kernel,Time,ResultType>::setThreshold( const podscalar thresh )
{
	__thresh = thresh;
}

template<class kernel,class Time,class ResultType>
typename GeneralSolver<kernel,Time,ResultType>::podscalar GeneralSolver<kernel,Time,ResultType>::threshold() const
{
	return __thresh;
}

template<class kernel,class Time,class ResultType>
void GeneralSolver<kernel,Time,ResultType>::setEstimatedError(const podscalar error )
{
	__estimated_error = error;
}

template<class kernel,class Time,class ResultType>
typename GeneralSolver<kernel,Time,ResultType>::podscalar GeneralSolver<kernel,Time,ResultType>::estimatedError() const
{
	return __estimated_error;
}

//////////////////End of impelementation of 'GeneralSolver'

};

#endif