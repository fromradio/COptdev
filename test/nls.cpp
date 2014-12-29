#include <Header>

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
	Function func(2,2);
	Vector initial_x(2);
	initial_x[0]=0;
	initial_x[1]=0;

	problem p(func, 2, 2);
	psolver lm(p, initial_x);
	lm.solve();
	lm.printInfo();
}