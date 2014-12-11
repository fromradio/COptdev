#include <Header>
 
typedef double		 					FT;
typedef COPT::Array<FT,int> 				Array;
typedef COPT::VectorBase<FT,int>			Vector;
typedef COPT::MatrixBase<FT,int>	 		Matrix;
typedef COPT::KernelTrait<FT,int>			kernel;
 

int main(int argc,char* argv[])
{
	Matrix A(4,2); 
 	A(0,0) = 1; A(0,1) = -1;
 	A(1,0) = -1;A(1,1) = 1; 
 	A(2,0) = 2; A(2,1) = -2;
 	A(3,0) = -3;A(3,1) = 1;
 	Vector b(4);
 	b[0] = 1;
 	b[1] = 2;
 	b[2] = 3;
 	b[3] = 4;
 	std::cout<<A<<std::endl;
 	std::cout<<b<<std::endl;
 	typedef COPT::LeastSquaresSolver<FT>       LeastSquares;
 	LeastSquares ls(A,b);
 	Vector x(2);
 	ls.solve(x);
 	ls.printInfo();
 	ls.setType(LeastSquares::LS);
 	ls.solve(x);
 	ls.printInfo();
 	ls.setType(LeastSquares::RLS);
 	ls.solve(x);
 	ls.printInfo();
}