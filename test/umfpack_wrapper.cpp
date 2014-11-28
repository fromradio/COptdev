#include <Header>

typedef double		 					FT;
typedef COPT::Array<FT> 				Array;
typedef COPT::VectorBase<FT,long>			Vector;
typedef COPT::MatrixBase<FT,long>	 		Matrix;
typedef COPT::SpMatrixBase<FT,long>			SpMatrix;
typedef COPT::UMFLinearSolver<SpMatrix>	UMFPackSolver;

int main(int argc , char* argv[])
{
	long rows = 10;
	long cols = 10;
	long elesize = 10;
	long* rowind = new long[elesize];
	long* colptr = new long[cols+1];
	FT*	vals = new FT[elesize];
	for ( int i = 0 ; i < elesize ; ++ i )
	{
		rowind[i]=i;
		vals[i]=2.0;
	}
	for ( int i = 0 ; i <= cols ; ++ i )
		colptr[i] = i;
	SpMatrix mat(rows,cols,elesize,colptr,rowind,vals);
	mat = 2*mat;

	UMFPackSolver solver(mat);
	// solver.printInfo();
	Vector v(10);
	v(0) = 1.0;
	std::cout<<solver.solve(v)<<std::endl;
	std::cout<<mat.solve(v)<<std::endl;
	delete[]rowind;
	delete[]colptr;
	delete[]vals;
}