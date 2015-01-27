#include "Core"

typedef COPT::MatrixBase<double,int> Matrix;

int main(int argc,char *argv[])
{
	Matrix mat = Matrix::random(5,2);
	std::cout<<mat<<std::endl;
	COPT::SVD<Matrix> svd(mat);
	std::cout<<mat<<std::endl; 
	std::cout<<svd.U()*svd.S()*svd.VT()<<std::endl;
}