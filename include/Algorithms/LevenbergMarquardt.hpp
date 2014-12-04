#ifndef LEVENBERGMARQUARDT_H
#define LEVENBERGMARQUARDT_H

#include <iostream>
#include <cmath>

namespace COPT
{
typedef VectorBase<double>			VT;
typedef MatrixBase<double>	 		Matrix;

VT LevenbergMarquardt(
	/*
	 *				the list of vector functions that are used
	 *				
	 */
	VectorFunctionSystem<VT> *vfs,
	//initial value
	VT initial_x,
	//function dimension
	int m,
	//variable dimension
	int n,
	//constant
	double tao=1e-3
	)
{
	
	/*
	 *	non-linear square solver
	 */

	/*
		Levenberg Marquardt method
	*/
	VT xf = initial_x;
	//VT f = vfs->FunctionValue(xf);
	//Matrix J = vfs->JacobiFun(xf);

	Matrix A = vfs->JacobiFun(xf).transpose()*vfs->JacobiFun(xf);	
	VT g = vfs->JacobiFun(xf)*vfs->FunctionValue(xf);

	double max = fabs(A(0,0));
	for(int i=0; i<n; i++){
		if(fabs(A(i,i))>max)
			max=A(i,i);
	}
	double u = tao*max;
	int v = 2;

	/*
	 *	the iteration terminals if the error is less than a threshold
	 *		or the iteration number is larger than a threshold
	 */
	int iternum = 0;
	double e = 1e-8;
	
	bool found = ( g.dot(g) <= e );
	double p;

	VT x;
	
	while ( !found && iternum < 100 ){
		++ iternum;
		Matrix B(A);
		VT k(g);	
		for(int i=0; i<n; i++){
			B(i,i) = A(i,i)+u;
			k[i]=-g[i];
		}
			
		VT h = B.solve(k);

		if ( h.dot(h) < e*(e+xf.dot(xf)) )
		{
			found = true;
		}
		else
		{
			x = xf + h;
			//p = (F(xf)-F(x))/L(h,u,g);
			p = (vfs->FunctionValue(xf).dot(vfs->FunctionValue(xf))/2-vfs->FunctionValue(x).dot(vfs->FunctionValue(x))/2)/(0.5*h.dot(u*h-g));
			if (p > 0.0)
			{
				xf = x;
				A = vfs->JacobiFun(xf).transpose()*vfs->JacobiFun(xf);
				g = vfs->JacobiFun(xf).transpose()*vfs->FunctionValue(xf);
				found = ( g.dot(g) <= e );
				double uu=1/3;
				if ( 1-(2*p-1)*(2*p-1)*(2*p-1) > uu )
				{
					uu = 1-(2*p-1)*(2*p-1)*(2*p-1);
				}
				u = u*uu;
				v = 2;
			}
			else
			{
				u = u*v;
				v = 2*v;
			}
					
		}
	}
	std::cout<<"The iteration number is "<<iternum<<std::endl;
	return xf;
}

} // End of namespace COPT

#endif