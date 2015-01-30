#include "Core"

typedef COPT::KernelTrait<double> 		kernel;
typedef kernel::SpMatrix 				SpMatrix;
typedef SpMatrix::Triplet 				Triplet;

typedef kernel::Matrix 				Matrix;
typedef kernel::Vector 				Vector;
typedef COPT::APGSolver<kernel, COPT::SolverTimeStatistics> APGSolver;


const int m = 7;
const int r = 1;
const int n = 10;
const int s = 9; // random element size

int main(int argc, char* argv[])
{
	// kernel::Matrix U,V;
	// U = kernel::Matrix::random(m,r);
	// V = kernel::Matrix::random(r,n);
	// kernel::Matrix T = U*V;
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
  

	Matrix A1,A2,A,D;
	int lam = 0.5;
	A1 = Matrix::random(m,r);
	A2 = Matrix::random(r,n);
	A = A1*A2*100;


	std::vector<int> aa(m*n);
	for ( int i = 0 ; i < aa.size() ; ++i ) aa[i]=i;
	std::random_shuffle(aa.begin(),aa.end());
	kernel::SpMatrix spE;
	std::vector<Triplet> tris;
	for ( int i = 0 ; i < s ; ++ i )
	{
		tris.push_back(Triplet(aa[i]-(aa[i]/m)*m,aa[i]/m,A(aa[i]-(aa[i]/m)*m,aa[i]/m)));
	}
	spE.setFromTriplets(m,n,tris.begin(),tris.end());
	Matrix E = spE.toDenseMatrix();
	std::cout<<"the fucking E?"<<std::endl;
	std::cout<<E<<std::endl;  
 
    // Matrix  E = Matrix::random(m,n);

    // E.resize(m,n);
    // E(0,0) = 0.01;
    // E(1,0) = 0.02;
    // E(2,2) = 0.03;
    // E(3,5) = 0.1;
    // E(4,8) = 2;
    // E(9,3) = 1;
	// std::cout<<"E:"<<std::endl;
	// std::cout<<E<<std::endl;

    D = A + E;

    std::cout<<"the shit D:"<<std::endl;
    std::cout<<D<<std::endl;
 
	APGSolver apg(D,lam);
	apg.solve();

	std::cout<<"the original A:"<<std::endl;
	std::cout<<A<<std::endl;

	// std::cout<<"error"<<std::endl;
	// std::cout<<(A - apg.result())<<std::endl;
	// std::cout<<"resultA:"<<std::endl;
	// std::cout<<apg.result()<<std::endl;
	// std::cout<<"resultE:"<<std::endl;
	// std::cout<<apg.result_E()<<std::endl;
}

