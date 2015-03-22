#include "Core"
#include <Eigen/Dense>
#include <iostream>
#include <time.h>

typedef Eigen::MatrixXd Matrix;
typedef COPT::KernelTrait<double> kernel;
typedef kernel::Matrix 			CMatrix;

int main(int argc,char *argv[])
{
	Matrix m1(100,100),m2(100,100),m3(100,100);
	m1.setRandom();
	m2.setRandom();
	m3.setRandom();
	// std::cout<<m1<<std::endl;
	clock_t t1,t2;
	t1 = clock();
	for (int i = 0 ; i < 1000 ; ++ i )
	{
		m1 = m2*m3;
	}
	t2 = clock();
	std::cout<<t2-t1<<" clocks"<<std::endl;
	t1 = clock();
	for (int i = 0 ; i < 1000 ; ++ i )
	{
		m2*m3;
	}
	t2 = clock();
	std::cout<<t2-t1<<" clocks"<<std::endl;
	CMatrix mm1(100,100),mm2(100,100),mm3(100,100);
	mm1 = CMatrix::random(100,100);
	mm2 = CMatrix::random(100,100);
	mm3 = CMatrix::random(100,100);

	CMatrix mm2t = mm2.transpose();
	t1 = clock();
	for (int i = 0 ; i < 1000 ; ++ i )
	{
		mm1 = mm2*mm3;
	}
	t2 = clock();
	std::cout<<t2-t1<<" clocks"<<std::endl;
	t1 = clock();
	for (int i = 0 ; i < 1000 ; ++ i )
	{
		m1 = m2.transpose()*m2;
	}
	t2 = clock();
	std::cout<<t2-t1<<" clocks"<<std::endl;
	for (int i = 0 ; i < 1000 ; ++ i )
	{
		mm1 = mm2.transpose()*mm2;
	}
	t2 = clock();
	std::cout<<t2-t1<<" clocks"<<std::endl;
	t1 = clock();
	for (int i = 0 ; i < 1000 ; ++ i )
	{
		m1 = m2.transpose()*m3;
	}
	t2 = clock();
	std::cout<<t2-t1<<" clocks"<<std::endl;
	for (int i = 0 ; i < 1000 ; ++ i )
	{
		mm2t*mm3;
	}
	t2 = clock();
	std::cout<<t2-t1<<" clocks"<<std::endl;
}