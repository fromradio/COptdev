#include "Core"

typedef COPT::KernelTrait<double,int>				kernel;

typedef kernel::Matrix	 		Matrix;
int main(int argc, char* argv[])
{
	kernel::Matrix m=kernel::Matrix::random(3000,4000);
	clock_t c1,c2;
	c1 = clock();
	m.operationNorm();
	c2 = clock();
	std::cout<<(double)(c2-c1)/CLOCKS_PER_SEC<<std::endl;
}