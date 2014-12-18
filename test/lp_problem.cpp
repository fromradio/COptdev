#include <Header>

typedef COPT::KernelTrait<double,int> 		kernel;
typedef kernel::Vector 						Vector;
typedef kernel::Matrix 						Matrix;
typedef COPT::LinearConstraint<kernel>		constraint;
typedef COPT::BPProblem<kernel> 			problem;
typedef COPT::BPSolver<problem,COPT::SolverTimeStatistics>				solver;

const int m = 2;
const int n = 5;

int main(int argc,char *argv[])
{
	Matrix A = Matrix::random(m,n);
	Vector b = Vector::random(m);

	/** x + y = 1.0, min |x|+|y| */
	// Matrix A(1,2);
	// A(0,0) = 1.0;
	// A(0,1) = 1.0;
	// Vector b(1); 
	// b(0) = 1;
	// std::cout<<"A is "<<std::endl<<A<<std::endl;

	problem p(n,A,b);
	p.validation();
	std::cout<<"p matrix is "<<std::endl<<p.matA()<<std::endl;
	solver sol(p);
	sol.solve();

	// Vector v(n);
	// v(0) = 0.1;
	// v(1) = -1.0;
	// std::cout<<p.objective(v)<<std::endl;
	// constraint con(2);
	// Vector vec(2);
	// vec(0)=1;
	// std::cout<<con.feasible(vec)<<std::endl;
}