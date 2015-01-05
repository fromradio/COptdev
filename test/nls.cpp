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

//	Function func(2,2);
//	Vector initial_x(2);
//	initial_x[0]=1;
//	initial_x[1]=2;

	Function func(4,4);
	Vector initial_x(4);
	initial_x[0]=0;
	initial_x[1]= -1.5;
	initial_x[2]= -1;
	initial_x[3]=1.5;

//	problem p(func, 2, 2);
	problem p(func, 4, 4);
	psolver lm(p, initial_x);
	lm.solve();
	lm.printInfo();
	std::cout<<"The objective value is  "<<lm.objective()<<std::endl;
}