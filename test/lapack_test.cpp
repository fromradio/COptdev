#include "Core"
#define size 3

typedef COPT::KernelTrait<double> 		kernel;
typedef kernel::Vector 					Vector;
typedef kernel::Matrix 					Matrix;

int main(int argc,char* argv[])
{
	// int i,j,c1,c2,pivot[size],ok;
	// float A[size][size],b[size],AT[size*size];
 
	Matrix m(3,3);

	m(0,0) = 1.0; m(0,1) = 1; m(0,2) = -1.0;
	m(1,0) = 0.0; m(1,1) = 1.0; m(1,2) = 1.0;
	m(2,0) = 1.0; m(2,1) = 0.0; m(2,2) = -1.0;

	Vector b(size);
	b[0] = -1.3;
	b[1] = -0.1;
	b[2] = 1.8;
	std::cout<<m<<std::endl;
	std::cout<<m.lapackSolve(b)<<std::endl;
	std::cout<<"b is "<<b<<std::endl;
	std::cout<<m.solve(b)<<std::endl;
	std::cout<<m.leastSquareSolve(b)<<std::endl;

	Vector v = m.solve(b);
	Vector  vv = m.lapackSolve(b);
	std::cout<<m*v<<std::endl;
	std::cout<<"vv is "<<vv<<std::endl;
	std::cout<<m*vv<<std::endl;
	std::cout<<b<<std::endl;
}