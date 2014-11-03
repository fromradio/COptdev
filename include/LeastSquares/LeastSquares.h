#ifndef LEASTSQUARE_H
#define LEASTSQUARE_H

// Implemented by Ruimin Wang
#include "../Frame/OptFrame.h"

namespace COpt
{
class LeastSquareSolver:public VectorOptimization
{
public:

	// type 
	enum LeastSquareType
	{
		Tra,
		LMS,
		RLS
	};
	// Constructor
	LeastSquareSolver();
	LeastSquareSolver(int dim);

	// setter and getter
	void setDimension(int dim);
	// solve least squares
	void solve(const Matrix& A,const Vector& b);

	/*
		Traditional Least Squares
	*/
	// solve least square problem using traditional method
	void solveTraditional(const Matrix&A,const Vector& b);



	/*
		Least Mean Squares method
	*/
	// solve least sqaures using least mean squares
	void solveLMS(const Matrix& A,const Vector& b , double length = 0.01);
	// initialize least mean squares method
	void initializeLMS();
	void initializeLMS( int dim );
	void initializeLMS( int dim , const Vector& v);
	// update least mean squares by single element
	void updateLMS(const Vector& a,double d);



	/*
		RLS method
	*/
	// solve least squares using RLS
	void solveRLS(const Matrix& A,const Vector& b);


private:
	// the type of the solver
	LeastSquareType __lstype;
	// the dimension of the solver
	int __dim;
	// step length of least mean squares
	double __lsm_length;
};
};

#endif