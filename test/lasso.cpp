#include <Header>    
typedef double 		FT;
typedef COPT::KernelTrait<FT>		kernel;
typedef kernel::Matrix 				Matrix;
typedef kernel::Vector 				Vector;
typedef COPT::LassoProblem<kernel>	problem;
typedef COPT::LassoProximalSolver<problem>	solver;

const int m = 5000;
const int n = 2000;
int main(int argc,char*argv[])
{
	Matrix A = Matrix::random(m,n);
	Vector b = Vector::random(m);
	///b(3) = 2.0;
	//b(4) = 1.0;
	problem sol(A,b,0.2);
	solver pro(sol,0.5,0.7,100);
	pro.solve();
}