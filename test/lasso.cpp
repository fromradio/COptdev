#include <Header>    
#include <IO>
typedef double 		FT;
typedef COPT::KernelTrait<FT>		kernel;
typedef kernel::Matrix 				Matrix;
typedef kernel::Vector 				Vector;
typedef COPT::LassoProblem<kernel>	problem;
typedef COPT::LassoADMMSolver<problem,COPT::SolverTimeStatistics>	solver;
typedef COPT::LassoProximalSolver<problem>	psolver;

const int m = 50;
const int n = 100;
int main(int argc,char*argv[])
{
	Matrix A;
	Vector b;
	A = Matrix::random(m,n);
	// A = Matrix(2,2);
	// A(0,0) = 1.0; A(1,1)=1.0;
	b = Vector::random(m);
	// b = Vector(2);
	// b(0)=1.0;b(1)=1.0;
	///b(3) = 2.0;
	//b(4) = 1.0;
	// COPT::readMtxFile("data/A1.mtx",A);
	// COPT::readMtxFile("data/b1.mtx",b);
	clock_t t1,t2;
	t1 = clock();
	problem pro(A,b,0.4);
	solver sol(pro,0.5,100000);
	

	
	t2 = clock();
	sol.solve();
	sol.printTimeInfo();
	psolver psol(pro,0.5,1.0,1000000);
	psol.solve();
}