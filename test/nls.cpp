#include "Core"

typedef double		 				FT;
typedef COPT::Array<FT,int> 				Array;
typedef COPT::VectorBase<FT,int>			Vector;
typedef COPT::MatrixBase<FT,int>	 		Matrix;
typedef COPT::KernelTrait<FT,int>			kernel;

typedef COPT::NonLinearSquareProblem<kernel>	problem;
typedef COPT::NonLinearSquare<problem,COPT::SolverTimeStatistics>	psolver;

int main(int argc,char* argv[])
{ 
	typedef COPT::VectorFunctionSystem<FT>	Function;
	Function func(10,2);
	Vector initial_x(2);
	initial_x[0]=0.3;
	initial_x[1]=0.4;
	problem p(func);
	psolver lm(p, initial_x);
	lm.solve();
	std::cout<<"The dimension is (m, n) = ("<<p.functionDimension()<<". "<<p.variableDimension()<<")"<<std::endl;
	lm.printInfo();
	std::cout<<"The objective value is  "<<lm.objective()<<std::endl;
}