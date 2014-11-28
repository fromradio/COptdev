#include <Header>
#include <IO>

typedef double		 					FT;
typedef COPT::Array<FT> 				Array;
typedef COPT::VectorBase<FT>			Vector;
typedef COPT::MatrixBase<FT>	 		Matrix;
typedef COPT::SpMatrixBase<FT>			SpMatrix;
typedef COPT::UMFLinearSolver<SpMatrix>	UMFPackSolver;

void initrand()
{
    srand((unsigned)(time(0)));
}

double randdouble()
{
    return rand()/(double(RAND_MAX)+1);
}

int main(int argc , char* argv[])
{
	std::cout<<argc<<std::endl;
	if( argc == 0 )
		return 0;
	else if (argc == 2)
	{
		SpMatrix m;
		Vector vec;
		readMtxFile(argv[1],m);
		std::cout<<"m's element size is "<<m.elementSize()<<std::endl;
		Vector v;
		v.resize(m.rows());
		initrand();
		for ( int i = 0 ; i < m.rows() ; ++ i )
			v[i] = randdouble();
		vec = m*v;
		std::cout<<"v is "<<v<<std::endl;
		Vector result = m.solve(vec);
		std::cout<<"result is "<<result<<std::endl;
		std::cout<<"error is "<<std::sqrt((result-v).squaredNorm())<<std::endl;
	}
	else if ( argc == 3 )
	{
		SpMatrix m;
		Vector vec;
		readMtxFile(argv[2],vec); 
		readMtxFile(argv[1],m);
		std::cout<<"m's element size is "<<m.elementSize()<<std::endl;
		Vector result = m.solve(vec);
		std::cout<<"result is "<<result<<std::endl;
		std::cout<<"error is "<<std::sqrt((m*result-vec).squaredNorm())<<std::endl;
	}
	else
		return 0;
}