#include <Header>

typedef double		 					FT;
typedef COPT::Array<FT> 				Array;
typedef COPT::VectorBase<FT,long>		Vector;
typedef COPT::MatrixBase<FT,long>	 	Matrix;
typedef COPT::TripletBase<FT>		 	Triplet;
typedef COPT::SpMatrixBase<FT>			SpMatrix;
typedef COPT::KernelTrait<FT,long>		kernel;
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
	OMPSolver solver(A,A.transpose()*A);
	solver.solve(b,3 );
	std::cout<<"result is "<<solver.result()<<std::endl;  
	// OMPSolver::normalizeAtomMatrix(A,norms);
	// std::cout<<A<<std::endl;
	// std::cout<<"norms are "<<norms<<std::endl;
	// Matrix ATA = A.transpose()*A;
	// OMPSolver::normalizeAtomMatrix(A,ATA,norms);
	// std::cout<<ATA<<std::endl;
	// std::list<long> indices;
	// indices.push_back(OMPSolver::findIndex(A,b));
	// Vector co;
	// Vector residual;
	// OMPSolver::updateCoefficient(A,A.transpose(),ATA,b,indices,co,residual);
	// std::cout<<co<<std::endl;
	// std::cout<<residual<<std::endl;
}