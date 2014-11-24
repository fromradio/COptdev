#include <Header>


typedef double		 					FT;
typedef COPT::Array<FT> 				Array;
typedef COPT::VectorBase<FT>			Vector;
typedef COPT::MatrixBase<FT>	 		Matrix;
typedef COPT::KernelTrait<FT>			kernel;
typedef COPT::SpMatrixBase<FT>			SpMatrix;

int main(int argc,char* argv[])
{
	size_t rows = 10;
	size_t cols = 10;
	size_t elesize = 10;
	size_t* rowind = new size_t[elesize];
	size_t* colptr = new size_t[cols+1];
	FT*	vals = new FT[elesize];
	for ( int i = 0 ; i < elesize ; ++ i )
	{
		rowind[i]=i;
		vals[i]=1.0;
	}
	for ( int i = 0 ; i <= cols ; ++ i )
		colptr[i] = i;
	SpMatrix mat(rows,cols,elesize,rowind,colptr,vals);
	std::cout<<mat(0,1)<<std::endl;
	std::cout<<mat(1,1)<<std::endl;  
	std::cout<<mat(9,9)<<std::endl;
}