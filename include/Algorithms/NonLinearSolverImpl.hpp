//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

namespace COPT
{
template<class VFunc>
NonLinearSolver<VFunc>::NonLinearSolver(
	const VFunc& func,
	const typename VFunc::ScalarType tol,
	const typename VFunc::ScalarType maxiternum,
	const typename VFunc::ScalarType c1,
	const typename VFunc::ScalarType c2,
	const typename VFunc::ScalarType ratio,
	const typename VFunc::ScalarType sigma,
	const int tracknum,
	const typename NonLinearSolver<VFunc>::SolverType type
	)
	:
	__func(func),
	__tol(tol),
	__maxiternum(maxiternum),
	__c1(c1),
	__c2(c2),
	__ratio(ratio),
	__sigma(sigma),
	__tracknum(tracknum),
	__status(false),
	__type(type)
{
}

/*		set the type of the solver: SDM, NM or BFGS
 *
 */
template<class VFunc>
void NonLinearSolver<VFunc>::setType(const SolverType type)
{
	__type = type;
}

/*		set the maximum iteration number
 *		/param maxiternum: 			the maximum itration number
 */
template<class VFunc>
void NonLinearSolver<VFunc>::setIterationNum(
	const int maxiternum)
{
	__maxiternum = maxiternum;
}

/* 		set the error threshold
 *		/param tol:					the error threshold
 */
template<class VFunc>
void NonLinearSolver<VFunc>::setErrorThreshold(
	const typename VFunc::ScalarType tol
	)
{
	__tol = tol;
}

/*		Kernel function
 *		Solving the problem
 */
template<class VFunc>
bool NonLinearSolver<VFunc>::doSolve()
{
	switch(__type)
	{
	case SDM:
	{
		steepestDescentUsingBackTracking(
			__func,
			__ratio,
			__c1,
			__x,
			__error,
			__iters,
			__tracknum);
	}
	break;
	case NM:
	{
		newtonMethod(
			__func,
			__x,
			__error,
			__iters);
	}
	break;
	case BFGS:
	{
		BFGSMethod(
			__func,
			__c1,
			__c2,
			__sigma,
			__x,
			__error,
			__iters,
			__ratio,
			__tracknum);
	}
	break;
	default:
	{
		throw COException("Unkown type of Nonlinear Solver! Please check the code!");
	}
	break;
	}
	if(__error<__tol){
		return true;
	}
	else{
		return false;
	}
}

/*		Solver interface for users
 *		/param x:			inital point
 */
template<class VFunc>
void NonLinearSolver<VFunc>::solve(
	const typename VFunc::Vector& x)
{
	__x = x;
	__error = __tol;
	__iters = __maxiternum;
	__status = doSolve();
}

template<class VFunc>
void NonLinearSolver<VFunc>::printInfo()
{
	std::cout<<"The result of non-linear solver by ";
	switch(__type){
	case SDM:
	{
		std::cout<<"Steepest descent method";
	}
	break;
	case NM:
	{
		std::cout<<"Newton's method";
	}
	break;
	case BFGS:
	{
		std::cout<<"BFGS method";
	}
	break;
	default:
	{
		throw COException("Unkown type of non-linear solver");
	}
	break;
	}
	std::cout<<std::endl;
	if(__status){
		std::cout<<"An at least local minimal has been successfully found as: "<<std::endl<<__x<<std::endl;
		std::cout<<"The whole procedure takes "<<__iters<<" iterations with norm of "<<__error<<std::endl;
	}
	else
	{
		std::cout<<"Solver fails to find minimal in "<<__iters<<" step. The norm of gradient of the function currently is "<<__error<<std::endl;
	}
}

} // End of namespace COPT