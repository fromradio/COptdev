//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#include <Header>


typedef COPT::KernelTrait<double,int> 			kernel;
typedef kernel::Vector 							Vector;
typedef kernel::Matrix 							Matrix;
typedef COPT::LU<Matrix> 						LU;

const int m = 5;
const int n = 2;
int main( int argc , char *argv[] )
{
	Matrix mat = Matrix::random(m,m);
	Eigen::MatrixXd mm(m,m);
	for (int i = 0 ; i < m ; ++ i )
		for ( int j = 0 ; j < m ; ++ j )
			mm(i,j) = mat(i,j);
	Vector v = Vector::random(m);
	Eigen::VectorXd vec(m);
	for ( int i = 0 ; i < m ; ++ i )
		vec(i) = v(i);
	Matrix rhb = Matrix::random(m,n);
	// std::cout<<"matrix is "<<std::endl<<mat<<std::endl;
	LU lu(mat);
	std::cout<<lu.solve(rhb)<<std::endl;
	std::cout<<lu.inverse()<<std::endl;
	std::cout<<"result of eigen"<<std::endl<<mm.inverse()<<std::endl;
	// std::cout<<Eigen::PartialPivLU<Eigen::MatrixXd>(mm).solve(vec)<<std::endl;
}