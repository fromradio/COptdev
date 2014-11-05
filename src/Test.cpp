#include <iostream>
#include <cmath>
#include "Function.h"
#include "Vector.h"
#include "Matrix.h"
#include "NumDifferential.h"
#include "Solver.h"
#include "FuncRepository.h"
#include "Conjugate.h"

typedef double 			FT;
typedef COPT::Vector<FT>			Vector;
typedef COPT::Matrix<FT>			Matrix;

int main(int argc,char* argv[])
{
	/*
	 * conjugate gradient method test
	 */
	Vector b(3),x(3);
	Matrix A(3,3);
	A(0,0) = 1;
	A(1,1) = 1;
	A(2,2) = 1;

	b(0) = 1;b(1) = 1;b(2) = 1;

	int iters = 100;
	FT tol = 1e-6;
	std::cout<<x<<std::endl;
	COPT::conjugateGradientWithoutPrecondition(A,b,x,iters,tol);
	std::cout<<"iteration used: "<<iters<<std::endl;
	std::cout<<"final error estimated "<<tol<<std::endl;
	std::cout<<x<<std::endl;
	// double* data = new double[10];
	// for ( int i = 0 ; i < 10 ; ++ i )
	// 	data[i] =  1.0;
	// std::cout<<data[0]<<' '<<"1"<<std::endl;
	// COPT::Vector<double> vec(10);
	// COPT::Vector<double> v(vec);
	// std::cout<<vec.size()<<std::endl;
	// // vec[0] = 1.0;
	// std::cout<<"1"<<std::endl;
	// vec[0] = 2.0;
	// v[1] = 1.0;
	// COPT::Vector<double> vv = vec+v;
	// std::cout<<vv[0]<<std::endl;
	// std::cout<<vv(1)<<std::endl;
	// const COPT::Vector<double> vvv(vv);
	// std::cout<<vvv(2)<<std::endl;
	// COPT::Matrix<double> mat(2,3);
	// mat(0,1) = 2;
	// std::cout<<mat(0,1)<<std::endl;
	// std::cout<<vec[0]<<' '<<v[-1]<<" Hello World!"<<std::endl;
	// COPT::QuadFunction<double> cf(1.0,-2.0); 
	// COPT::RootSolver<COPT::ScalarFunction<double> > rs(cf);
	// std::cout<<rs.solve();
	// typedef COPT::Vector<double> Vector;

	// COPT::Vector<double> vec(2);
	// vec[0] = 1.0;
	// vec[1] = 2.0;
	// COPT::VectorCosineFunction<Vector> cf(vec);
	// Vector v(2);
	// v[0] = 1.57;
	// v[1] = 0;
	// std::cout<<cf(v)<<std::endl;
	// std::cout<<cf.gradient(v)<<std::endl; 
	// std::cout<<rs.
	// std::cout<<cf(1.57)<<std::endl;
	// COPT::ScalarDifferential<COPT::CosineFunction<double> > diff(cf);	
	// std::cout<<cf.diff(1.57)<<' '<<diff.diff(1.57)<<' '<<fabs(cf.diff(1.57)-diff.diff(1.57))<<std::endl;
	// COPT::SAFE_DELETE(data);
}