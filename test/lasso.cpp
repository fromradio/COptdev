#include <Header>    
typedef double 		FT;
typedef COPT::KernelTrait<FT>		kernel;
typedef kernel::Matrix 				Matrix;
typedef kernel::Vector 				Vector;
typedef COPT::LassoProblem<kernel>	problem;
typedef COPT::LassoADMMSolver<problem>	solver;
typedef COPT::LassoProximalSolver<problem>	psolver;

const int m = 500;
const int n = 200;
int main(int argc,char*argv[])
{
	Matrix A = Matrix::random(m,n);
	Vector b = Vector::random(m);
	///b(3) = 2.0;
	//b(4) = 1.0;
	problem pro(A,b,0.2);
	solver sol(pro,0.5,1000);
	psolver psol(pro,0.5,0.7,1000);
	sol.solve();
	psol.solve();
}