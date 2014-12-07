#include <Header>    
typedef double 		FT;
typedef COPT::KernelTrait<FT>		kernel;
typedef kernel::Matrix 				Matrix;
typedef kernel::Vector 				Vector;
typedef COPT::LassoProblem<kernel>	problem;
typedef COPT::LassoProximalSolver<problem>	solver;

int main(int argc,char*argv[])
{
	Matrix A = Matrix::random(5000,2000);
	Vector b = Vector::random(5000);
	///b(3) = 2.0;
	//b(4) = 1.0;
	problem sol(A,b,0.2);
	solver pro(sol);
	pro.solve();
	//  
}