#include "Core"

typedef COPT::KernelTrait<double> 		kernel;
typedef kernel::SpMatrix 				SpMatrix;
typedef SpMatrix::Triplet 				Triplet;


typedef kernel::Matrix 				Matrix;
typedef kernel::Vector 				Vector;
typedef COPT::APGSolver<kernel, COPT::SolverTimeStatistics> APGSolver;


const int m = 5;
const int r = 1;
const int n = 8;
const int s = 7; // random element size
 
int main(int argc, char* argv[])
{
	// kernel::Matrix U,V;     
	// U = kernel::Matrix::random(m,r);
	// V = kernel::Matrix::random(r,n);
	// kernel::Matrix T = U*V; 
	// addSparseNoise(T,s,0.01);
	// std::cout<<T<<std::endl;
	// COPT::LowRankALM<kernel::Matrix> alm(T,0.5);
	// alm.solve();
	// std::cout<<"A is "<<std::endl<<alm.A()<<std::endl;
	// std::cout<<"E is "<<std::endl<<alm.E()<<std::endl;

	// std::vector<int> tt(m*n);
	// for ( int i = 0 ; i < tt.size() ; ++i ) tt[i]=i; 
	// std::random_shuffle(tt.begin(),tt.end());
	// kernel::SpMatrix spmat;
	// std::vector<Triplet> tris;
	// for ( int i = 0 ; i < s ; ++ i )
	// {
	// 	tris.push_back(Triplet(tt[i]-(tt[i]/m)*m,tt[i]/m,T(tt[i]-(tt[i]/m)*m,tt[i]/m)));
	// }
	// spmat.setFromTriplets(m,n,tris.begin(),tris.end());
	// COPT::APGLSolver<kernel> solver(spmat);
	// solver.solve();
	// std::cout<<spmat.toDenseMatrix()<<std::endl;
   

	Matrix U,V;
	U = kernel::Matrix::random(m,r);
	V = kernel::Matrix::random(r,n);
	Matrix A = U*V; 
	Matrix D = A;
	addSparseNoise(D,s,0.1);
	Matrix E = D - A;
	APGSolver apg(D,0.5); 
	apg.solve();
	std::cout<<"the fucking result E:"<<std::endl;
	std::cout<<apg.result_E()<<std::endl;
	std::cout<<"oh hehe this is the real E:"<<std::endl;
	std::cout<<E<<std::endl;
}


