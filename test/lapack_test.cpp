#include <Header>
#define size 3

typedef COPT::KernelTrait<float> 		kernel;
typedef kernel::Vector 					Vector;
typedef kernel::Matrix 					Matrix;

int main(int argc,char* argv[])
{
	int i,j,c1,c2,pivot[size],ok;
	// float A[size][size],b[size],AT[size*size];

	Matrix m(3,3);

	m(0,0) = 1.0; m(0,1) = 1; m(0,2) = -1.0;
	m(1,0) = 0.0; m(1,1) = 1.0; m(1,2) = 1.0;
	m(2,0) = 1.0; m(2,1) = 0.0; m(2,2) = -1.0;

	Vector b(size);
	b[0] = -1.3;
	b[1] = -0.1;
	b[2] = 1.8;

	// for ( i = 0 ; i < size ; ++ i )
	// {
	// 	for ( j = 0 ; j < size ; ++ j )
	// 		AT[j+size*i]=A[j][i];
	// }

	// c1 = size;
	// c2 = 1;

	// COPT::copt_lapack_gesv(&c1,&c2,AT,&c1,pivot,b,&c1,&ok);

	// Vector vec(c1,b);
	// Matrix mat(size,size,AT);
	// mat = Matrix::identity(size,size);

	std::cout<<m<<std::endl;

	clock_t start, end;
	start = clock();
	std::cout<<m.lapackSolve(b)<<std::endl;
	std::cout<<"b is "<<b<<std::endl;
	std::cout<<m.solve(b)<<std::endl;

	Vector v = m.solve(b);
	Vector  vv = m.lapackSolve(b);
	std::cout<<m*v<<std::endl;
	std::cout<<"vv is "<<vv<<std::endl;
	std::cout<<m*vv<<std::endl;
	std::cout<<b<<std::endl;

	for ( j = 0 ; j < size ; ++ j ) printf("%e\n",b[j]);
}