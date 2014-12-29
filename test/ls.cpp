#include <Header>
#include <IO>
 
typedef double 		FT;
typedef COPT::KernelTrait<FT>		kernel;
typedef kernel::Matrix 				Matrix;
typedef kernel::Vector 				Vector;
typedef COPT::LeastSquaresProblem<kernel>	problem;
typedef COPT::LeastMeanSquareSolver<problem,COPT::SolverTimeStatistics>	lmssolver;
typedef COPT::LeastSquareSolver<problem,COPT::SolverTimeStatistics>	lssolver;
typedef COPT::RecursiveLeastSquareSolver<problem,COPT::SolverTimeStatistics> 	rlssolver;

const int m = 180;
const int n = 360;
int main(int argc,char* argv[])
{
	Matrix A;
	Vector b;
	A = Matrix::random(m,n);
	b = Vector::random(m);
	problem pro(A,b);
	lmssolver lmssol(pro);

	COPT::Printer::printType(A,std::cout);

	lmssol.solve();

}