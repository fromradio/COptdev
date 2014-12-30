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
 
  
const int m = 100;
const int n = 10;

int main(int argc,char* argv[])
{
	Matrix A;
	Vector b;
	A = Matrix::random(m,n);
	b = Vector::random(m);
	std::cout<<"A : "<<std::endl<<A<<std::endl;
	std::cout<<"b : "<<std::endl<<b<<std::endl;
	problem pro(A,b);
	lmssolver lmssol(pro);
	lssolver lssol(pro);
	rlssolver rlssol(pro);
    
	// COPT::Printer::printType(A,std::cout);
	// COPT::Printer::printType(b,std::cout);

    std::cout<<"results of three solvers:"<<std::endl;
	lmssol.solve();
	lssol.solve();
	rlssol.solve();

}