#include "Core"
#include "IO"

typedef COPT::KernelTrait<std::complex<double>,int>		kernel;
typedef kernel::Vector 						Vector;
typedef kernel::Matrix 						Matrix;

int main( int argc , char *argv[])
{
	Matrix mat(5,5);
	Vector vec;
	COPT::readMtxFile("data/b.mtx",vec);
	COPT::readMtxFile("data/R.mtx",mat);

	// std::cout<<mat<<std::endl;
	std::cout<<mat.transpose()*vec<<std::endl; 
}