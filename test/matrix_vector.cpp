#include "Core"

typedef double		 					FT;
typedef COPT::Array<FT,int> 				Array;
typedef COPT::VectorBase<FT,int>			Vector;
typedef COPT::MatrixBase<FT,int>	 		Matrix;
typedef COPT::KernelTrait<FT,int>			kernel;

int main(int argc,char* argv[])
{
	Matrix mat(2,5);
	mat(1,0) = 2;
	mat(1,1) = 4;
	mat(0,3) = 3;
	mat(0,4) = 1;
	std::cout<<mat<<std::endl;
	std::cout<<mat.col(0)<<std::endl;
	std::cout<<mat.col(1)<<std::endl;
	std::cout<<mat.row(0)<<std::endl;
	std::cout<<mat.row(1)<<std::endl;
	Vector v1(2),v2(5);
	v1(1) = 3;
	std::cout<<mat.col(0)+mat.col(1)<<std::endl;
	mat.col(0)(0) = 3.0;
	std::cout<<5.0*mat.col(0)<<std::endl;
	std::cout<<mat<<std::endl;
}