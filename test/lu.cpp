//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#include <Header>


typedef COPT::KernelTrait<std::complex<double>,int> 			kernel;
typedef kernel::Vector 							Vector;
typedef kernel::Matrix 							Matrix;
typedef COPT::LU<Matrix> 						LU;
typedef COPT::CholeskySolver<Matrix>			Cholesky;
typedef COPT::QR<Matrix> 						QR;

const int m = 5;
const int n = 2;
int main( int argc , char *argv[] )
{
	Matrix mat = Matrix::random(m,m);
	// Eigen::MatrixXd mm(m,m);
	// for (int i = 0 ; i < m ; ++ i )
	// 	for ( int j = 0 ; j < m ; ++ j )
	// 		mm(i,j) = mat(i,j);
	Vector v = Vector::random(m);
	// Eigen::VectorXd vec(m);
	// for ( int i = 0 ; i < m ; ++ i )
		// vec(i) = v(i);
	Matrix rhb = Matrix::random(m,n);
	// std::cout<<"matrix is "<<std::endl<<mat<<std::endl;
	LU lu(mat);
	Matrix mtm;
	mat.mtm(mtm);
	Cholesky cho(mtm);
	QR qr(mat);
	std::cout<<lu.solve(rhb)<<std::endl;
	lu.compute(mtm);
	std::cout<<lu.solve(mat.transpose()*rhb)<<std::endl;
	std::cout<<qr.solve(mat.transpose()*rhb)<<std::endl;
	// std::cout<<lu.inverse()<<std::endl;
	std::cout<<"solving result"<<std::endl;
	std::cout<<cho.solve(mat.transpose()*rhb)<<std::endl;

	// LU lu2(mtm);
	// std::cout<<"result of eigen"<<std::endl<<mm.inverse()<<std::endl;
	// std::cout<<lu2.inverse()<<std::endl;
	// std::cout<<lu2.inverse()*mtm<<std::endl;

	// std::cout<<"inverse is "<<std::endl<<cho.inverse()<<std::endl;
	// std::cout<<cho.inverse()*mtm<<std::endl; 
	// std::cout<<Eigen::PartialPivLU<Eigen::MatrixXd>(mm).solve(vec)<<std::endl;
}