#include "Header"
#include "IO"

typedef COPT::KernelTrait<double,int>		kernel;
typedef kernel::Vector 						Vector;
typedef kernel::Matrix 						Matrix;

int main( int argc , char *argv[])
{
	Matrix mat(5,5);
	// mat(0,0)=1.0;
	// mat(3,0)=1.0;
	// COPT::writeMtxFile("data/test_m.mtx",mat);
	COPT::readMtxFile("data/test_m.mtx",mat);
	std::cout<<mat<<std::endl;
}