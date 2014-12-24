#include <Header>    
#include <IO>
typedef double 		FT;
typedef COPT::KernelTrait<FT>		kernel;
typedef kernel::Matrix 				Matrix;
typedef kernel::Vector 				Vector;
typedef COPT::LassoProblem<kernel>	problem;
typedef COPT::LassoADMMSolver<problem,COPT::SolverTimeStatistics>	solver;
typedef COPT::LassoProximalSolver<problem,COPT::SolverTimeStatistics>	psolver;
typedef COPT::FISTALassoSolver<problem,COPT::SolverTimeStatistics> 	fsolver;

const int m = 180;
const int n = 360;
int main(int argc,char*argv[])
{
	Matrix A;
	Vector b;
	A = Matrix::random(m,n);
	// A = Matrix(2,2);
	// A(0,0) = 1.0; A(1,1)=1.0; A(1,0) = 1.0;
	b = Vector::random(m);
	// b = Vector(2);
	// b(0)=1.0;b(1)=1.0;
	// std::cout<<"A is "<<A<<std::endl; 
	// std::cout<<"A*b is "<<A*b<<std::endl;
	///b(3) = 2.0;
	//b(4) = 1.0;
	// COPT::readMtxFile("data/A1.mtx",A);
	// COPT::readMtxFile("data/b1.mtx",b);
	// clock_t t1,t2;
	// t1 = clock();
	problem pro(A,b,0.4);
	psolver sol(pro);
	

	COPT::Printer::printType(A,std::cout);
	
	// t2 = clock();
	sol.solve();
	// sol.printTimeInfo();
	// fsolver fsol(pro);
	// fsol.solve();
	// psolver psol(pro,0.5,1.0,1000000);
	// psol.solve();
}