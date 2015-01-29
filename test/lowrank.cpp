#include "Core"

typedef COPT::KernelTrait<double> 		kernel;
typedef kernel::SpMatrix 				SpMatrix;
typedef SpMatrix::Triplet 				Triplet;

const int m = 10;
const int r = 2;
const int n = 15;
const int s = 20; // random element size

int main(int argc, char* argv[])
{
	kernel::Matrix U,V;
	U = kernel::Matrix::random(m,r);
	V = kernel::Matrix::random(r,n);
	kernel::Matrix T = U*V;
	std::vector<int> tt(m*n);
	for ( int i = 0 ; i < tt.size() ; ++i ) tt[i]=i;
	std::random_shuffle(tt.begin(),tt.end());
	kernel::SpMatrix spmat;
	std::vector<Triplet> tris;
	for ( int i = 0 ; i < s ; ++ i )
	{
		tris.push_back(Triplet(tt[i]-(tt[i]/m)*m,tt[i]/m,T(tt[i]-(tt[i]/m)*m,tt[i]/m)));
	}
	spmat.setFromTriplets(m,n,tris.begin(),tris.end());
	COPT::APGLSolver<kernel> solver(spmat);
	solver.solve();
	std::cout<<spmat.toDenseMatrix()<<std::endl;
}