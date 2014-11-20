#include <Header>

typedef double 		FT;
typedef COPT::VectorBase<FT>	Vector;
typedef COPT::MatrixBase<FT> 	Matrix;
typedef COPT::KernelTrait<FT> 	kernel;
typedef COPT:: SimplexSolver<kernel>	Simplex;

int main(int argc,char* argv[])
{
	Vector c(4);
	c[0] = -3;
	c[1] = -2;
	c[2] = 0;
	c[3] = 0;
	Matrix A(2,4);
	A(0,0)=1;
	A(0,1)=1;
	A(0,2)=1;
	A(1,0)=2;
	A(1,1)=0.5;
	A(1,3)=1;
	Vector b(2);   
	b[0]=5;
	b[1]=8;
	std::vector<size_t> indb,indn;
	indb.push_back(2);
	indb.push_back(3);
	indn.push_back(0);
	indn.push_back (1); 
	std::cout<<"current status is "<<Simplex::oneSimplexStep(A,b,c,indb,indn);
	std::cout<<"current status is "<<Simplex::oneSimplexStep(A,b,c,indb,indn);
	std::cout<<"current status is "<<Simplex::oneSimplexStep(A,b,c,indb,indn);
}