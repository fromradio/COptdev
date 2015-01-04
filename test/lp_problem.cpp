#include "Core"
#include <IO>


typedef COPT::KernelTrait<std::complex<double>,int> 		kernel;
typedef kernel::Vector 						Vector;
typedef kernel::Matrix 						Matrix;
typedef COPT::LinearConstraint<kernel>		constraint;
typedef COPT::BPProblem<kernel> 			problem;
typedef COPT::BPSolver<problem,COPT::SolverTimeStatistics>				solver;

const int m = 2;
const int n = 5;

int main(int argc,char *argv[])
{
	Matrix R ;//= Matrix::random(m,n);
	Vector b ;//= Vector::random(m);

	COPT::readMtxFile("data/R.mtx",R);
	COPT::readMtxFile("data/b.mtx",b);

	Matrix A,iden;
	iden = Matrix::identity(R.rows(),R.rows());
	iden = (-0.8)*iden;
	A.combineAlongColumn(R,iden);

	/** x + y = 1.0, min |x|+|y| */
	// Matrix A(1,2);
	// A(0,0) = 1.0;
	// A(0,1) = 1.0;
	// Vector b(1); 
	// b(0) = 1;
	// std::cout<<"A is "<<std::endl<<A<<std::endl;

	std::cout<<A.rows()<<" "<<A.cols()<<std::endl;
	problem p(A.cols(),A,b);
	p.validation();
	// std::cout<<"p matrix is "<<std::endl<<p.matA()<<std::endl;
	solver sol(p,1.0,500,1e-12);
	sol.solve();

	COPT::writeMtxFile("data/result.mtx",sol.result());

	// Vector v(n);
	// v(0) = 0.1;
	// v(1) = -1.0;
	// std::cout<<p.objective(v)<<std::endl;
	// constraint con(2);
	// Vector vec(2);
	// vec(0)=1;
	// std::cout<<con.feasible(vec)<<std::endl;
}