#include <Header>

typedef COPT::KernelTrait<std::complex<double>,int> 			kernel;
typedef kernel::Matrix 							Matrix;
typedef kernel::Vector 							Vector;
typedef COPT::PartialEigenSolver<Matrix>		EigenSolver;
int main(int argc,char* argv[])
{
	Matrix m(2,2);
	m(0,0)=2;m(1,1)=1;
	Matrix mtm;
	m.mtm(mtm);
	mtm(0,0) = std::complex<double>(1.0,0);
	mtm(0,1) = std::complex<double>(0,-1);
	mtm(1,0) = std::complex<double>(0,1);
	std::cout<<mtm<<std::endl;
	EigenSolver solver(mtm);
	std::cout<<solver.computeLargestEigenvalue()<<std::endl;
	std::cout<<m.operationNorm()<<std::endl;
	m = Matrix::random(5,3);
	double e = m.operationNorm();
	for ( int i = 0 ; i < 100000 ; ++ i )
	{
		Vector v = Vector::random(3);
		double f1 = std::sqrt((m*v).squaredNorm());
		double f2 = std::sqrt(v.squaredNorm())*e;
		if(f1>f2)
			std::cout<<"error"<<std::endl;
		// std::cout<<f1<<' '<<f2<<std::endl;
	}
}