#ifndef NON_LINEAR_SQUARE_H
#define NON_LINEAR_SQUARE_H


namespace COPT
{

template<class T>
class NonLinearSquare
{
private:
	typedef VectorBase<T>					VT;
	typedef MatrixBase<T>	 			Matrix;
	/*
	 *		vector functions that are used
	 *				
	 */
	VectorFunctionSystem<T>*		__vfs;
	//	point value
	VT 						  		  __x;
	//	function dimension
	int 							  __m;
	//	variable dimension
	int 							  __n;
	//	constants
	T 							 	__tao;
	int 							  __v;
	int 						__iternum;


	/*
	 *		non-linear square solver
	 */
	void           			    doSolve();

public:
	/*
	 *		Constructor
	 */ 
	NonLinearSquare(
		VectorFunctionSystem<T>* vfs,
		int m,
		int n,
		T tao = 1e-3,
		int v = 2,
		int iternum = 100
		);
	
	/*
	 *		Destructor
	 */ 
	~NonLinearSquare();

	void solve(VT& initial_x);

	// 	output
	void printInfo();
	// 	optimal point
	VT& OptimalPoint() {return __x;}

};

template<class T>
NonLinearSquare<T>::NonLinearSquare(
	VectorFunctionSystem<T>* vfs,
	int m,
	int n,
	T tao,
	int v,
	int iternum
	)
    :
    __vfs(vfs),
    __m(m),
    __n(n),
    __tao(tao),
    __v(v),
    __iternum(iternum)
{

}

template<class T>
NonLinearSquare<T>::~NonLinearSquare()
{

}

template<class T>
void NonLinearSquare<T>::doSolve()
{
	/*
		Levenberg Marquardt method
	*/

	/*	A is approximately equal to Hessian	matrix 	*/	
	Matrix A = __vfs->JacobiFun(__x).transpose()*(__vfs->JacobiFun(__x));	
	/*	g is approximately equal to Gradient	*/
	VT g = __vfs->JacobiFun(__x)*(__vfs->FunctionValue(__x));

	double max = fabs(A(0,0));
	for(int i=0; i<__n; i++){
		if(fabs(A(i,i))>max)
			max=A(i,i);
	}
	double u = __tao*max;

	/*
	 *	the iteration terminals if the error is less than a threshold
	 *		or the iteration number is larger than a threshold
	 */
	int iternum = 0;
	/*	error	*/
	double e = 1e-8;
	
	bool found = ( g.dot(g) <= e );
	double p;

	VT x;
	
	while ( !found && iternum < __iternum ){
		++ iternum;
		Matrix B(A);
		VT k(g);	
		for(int i=0; i<__n; i++){
			B(i,i) = A(i,i)+u;
			k[i]=-g[i];
		}
			
		/*	Solve linear equations	*/	
		VT h = B.solve(k);

		if ( h.dot(h) < e*(e+__x.dot(__x)) )
		{
			found = true;
		}
		else
		{
			x = __x + h;

			/*	Judgment flag	p = (F(__x)-F(x))/L(h,u,g)	*/
			p = (__vfs->FunctionValue(__x).dot(__vfs->FunctionValue(__x))/2-__vfs->FunctionValue(x).dot(__vfs->FunctionValue(x))/2)/(0.5*h.dot(u*h-g));

			if (p > 0.0)
			{
				__x = x;
				A = __vfs->JacobiFun(__x).transpose()*(__vfs->JacobiFun(__x));
				g = __vfs->JacobiFun(__x).transpose()*(__vfs->FunctionValue(__x));
				found = ( g.dot(g) <= e );
				double uu=1/3;
				if ( 1-(2*p-1)*(2*p-1)*(2*p-1) > uu )
				{
					uu = 1-(2*p-1)*(2*p-1)*(2*p-1);
				}
				u = u*uu;
				__v = 2;
			}
			else
			{
				u = u*__v;
				__v = 2*__v;
			}			
		}
	}
}

template<class T>
void NonLinearSquare<T>::solve(VT& initial_x)
{
	__x = initial_x;
	doSolve();
}

template<class T>
void NonLinearSquare<T>::printInfo()
{
	std::cout<<"The optimal point is "<<__x<<std::endl;
}


} // End of namespace COPT

#endif
