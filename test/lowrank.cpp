#include "Core"

typedef COPT::KernelTrait<double> 		kernel;
typedef kernel::SpMatrix 				SpMatrix;
typedef SpMatrix::Triplet 				Triplet;

const int m = 500;
const int r = 20;
const int n = 700;
const int s = 500; // random element size

int main(int argc, char* argv[])
{
	kernel::Matrix U,V;
	U = kernel::Matrix::random(m,r);
	V = kernel::Matrix::random(r,n);
	kernel::Matrix T = U*V;
	addSparseNoise(T,s,0.01);
	std::cout<<T<<std::endl;
	COPT::LowRankALM<kernel::Matrix> alm(T,0.5);
	alm.solve();
	std::cout<<"A is "<<std::endl<<alm.A()<<std::endl;
	std::cout<<"E is "<<std::endl<<alm.E()<<std::endl;
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
}