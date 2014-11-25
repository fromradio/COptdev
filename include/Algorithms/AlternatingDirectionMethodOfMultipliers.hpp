//		Copyright (C) Jingyuan Hu, hujy@mathu.cn
//		Copyright (C) MathU

#ifndef ADMM_H
#define ADMM_H

//Alternating Direction Mathod Of Multipliers

namespace COPT{

/*
 *		judge how to change the term of vector
 *      /param rho:       the constant ratio
 *      /param a:         the term of the vector
 */    
template<class Scalar>
Scalar Sfunction(const Scalar rho,const Scalar a)
{
	if(a>1/rho)
		return a-1/rho;
	else if(a<-1/rho)
		return a+1/rho;
	else{
		return 0;
	} 
}

/*   
 *		Alternating Direction Mathod Of Multipliers      
 *      /param A:          the coefficient matrix
 *      /param b:          constant vector
 *      /param rho:        the constant ratio
 *      /param u:          the quotient of the Lagrange param and rho
 *      /param x:          the weight that satisfizes the equation
 *      /param z:          the inital z on input and optimized z on output
 *      /param number:     the most number of iteration 
 *      /param e:          the least error of iteration
 */    
template<class Scalar>
void LeastAbsoluteDeviationMethod(
	const MatrixBase<Scalar>& A,
	const VectorBase<Scalar>& b,
	const Scalar rho,
	VectorBase<Scalar>& u,
	VectorBase<Scalar>& x,
	VectorBase<Scalar>& z,
	int& number,
	Scalar& e)
{
	int count=number;
	int m=0;
	int n=0;
	number=0;
	double t;
	double error;
	MatrixBase<Scalar> At;
	MatrixBase<Scalar> Az;
	VectorBase<Scalar> Zstay;
	VectorBase<Scalar> c;

	t=e*e;
	At=A.transpose();
	Az=At*A;


	while(number<count){
		c=At*(b+z-u);
		x=Az.solve(c);

		Zstay=A*x-b+u;
		n=Zstay.size();
		m=0;
		while(m<n){
			Zstay[m]=Sfunction(rho,Zstay[m]);
			m=m+1;
		}
		error=(z-Zstay).squaredNorm()+(A*x-z-b).squaredNorm();
		z=Zstay;

		u=u+A*x-z-b;

		number=number+1;

		if(error<t)
			break;

	}
}
}

#endif