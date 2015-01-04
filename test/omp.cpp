#include "Core"

typedef double		 					FT;
typedef COPT::Array<FT,int> 			Array;
typedef COPT::VectorBase<FT,int>		Vector;
typedef COPT::MatrixBase<FT,int>	 	Matrix;
typedef COPT::TripletBase<FT,int>		 	Triplet;
typedef COPT::SpMatrixBase<FT,int>			SpMatrix;
typedef COPT::KernelTrait<FT,int>		kernel;
typedef COPT::OMPSolver<kernel>			OMPSolver;

int main(int argc,char*argv[])
{
	Matrix A=Matrix::identity(5,5);
	
	Vector b(5);
	b(1) = 1.0;
	b(0) = 1.0;
	b(2) = 3.0;
	// OMPSolver solver(A);
	std::cout<<b<<std::endl; 
	A(0,0) = 3.0;
	A(0,1) = 1.0;
	std::cout<<"A is "<<std::endl<<A<<std::endl;
	Vector norms;
	OMPSolver solver(A,A.transpose()*A );
	solver.solve(b,3 );
	std::cout<<"result is "<<solver.result()<<std::endl; 
	std::cout<<"fitting error is "<<solver.fittingError()<<std::endl;
}
