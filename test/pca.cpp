#include "Core"


typedef COPT::KernelTrait<double> 			kernel;
typedef kernel::Matrix 						Matrix;

int main(int argc, char *argv[])
{
	// Matrix mat = Matrix::identity(5,5);
	// mat(0,4)=1.0;
	// mat(0,3)=1.0;
	// mat.setSymmetricFlag(true);
	// COPT::EigenSolver<Matrix> ev(mat);
	Matrix mat = Matrix::random(5,7);
	COPT::PCA<Matrix> pca(mat,4);
	std::cout<<"result is "<<pca.solve(mat)<<std::endl;
	// std::cout<<ev.eigenValue()<<std::endl;
	// std::cout<<ev.eigenVector()<<std::endl;
	// std::cout<<ev.eigenVector()*Matrix::diag(5,5,ev.eigenValue())*ev.eigenVector().transpose()<<std::endl;
	// std::cout<<mat*ev.eigenVector().col(0)<<std::endl;
	// std::cout<<ev.eigenValue()(0)*ev.eigenVector().col(0)<<std::endl;
}