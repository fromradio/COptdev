#include "Core"

typedef COPT::KernelTrait<double> 				kernel;
typedef kernel::Matrix 							Matrix;
typedef COPT::Solver<double,Matrix> 			Solver;

double func(const Matrix& mat)
{
	return mat.frobeniusNorm();
}

int main(int argc, char *argv[])
{
	Solver sol(&func);
	Matrix iden = Matrix::identity(4,4);
	std::cout<<sol.objective(iden)<<std::endl;
}