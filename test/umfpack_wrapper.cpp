//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU 

#include "Core"
#include <IO>

typedef double		 						FT;
typedef COPT::KernelTrait<FT> 	 			kernel;
typedef COPT::Array<FT,int> 				Array;
typedef COPT::VectorBase<FT,int>			Vector;
typedef COPT::MatrixBase<FT,int>	 		Matrix;
typedef COPT::SpMatrixBase<FT,int >			SpMatrix;
typedef COPT::UMFLinearSolver<SpMatrix>		UMFPackSolver;

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
		COPT::readMtxFile(argv[1],m);
		// std::cout<<"m's element size is "<<m.elementSize()<<std::endl;
		Vector v;
		v.resize(m.rows());
		initrand();
		for ( int i = 0 ; i < m.rows() ; ++ i )
			v[i] = randdouble();
		vec = m*v;
		// std::cout<<"v is "<<v<<std::endl;
		UMFPackSolver solver(m);
		std::cout<<"solving starting...."<<std::endl;
		Vector result = solver.solve(vec);
		// // std::cout<<"result is "<<result<<std::endl;
		// std::cout<<m.toDenseMatrix()<<std::endl; 
		std::cout<<"the matrix is with size "<<m.rows()<<" "<<m.cols()<<std::endl;
		std::cout<<"compute time is "<<solver.computeTime()<<"s."<<std::endl;
		std::cout<<"solve time is "<<solver.solveTime()<<"s."<<std::endl;
		// solver.printInfo(); 
		std::cout<<"error is "<<std::sqrt((result-v).squaredNorm())<<std::endl;
	}
	else if ( argc == 3 )
	{
		SpMatrix m;
		Vector vec;
		COPT::readMtxFile(argv[2],vec); 
		COPT::readMtxFile(argv[1],m);
		// std::cout<<"m's element size is "<<m.elementSize()<<std::endl;
		UMFPackSolver solver(m);
		Vector result = solver.solve(vec); 
		std::cout<<"the matrix is with size "<<m.rows()<<" "<<m.cols()<<std::endl;
		std::cout<<"compute time is "<<solver.computeTime()<<"s."<<std::endl;
		std::cout<<"solve time is "<<solver.solveTime()<<"s."<<std::endl;
		std::cout<<"error is "<<std::sqrt((m*result-vec).squaredNorm())<<std::endl;
	}
	else
		return 0;
}