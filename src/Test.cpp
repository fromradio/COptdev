#include <Header>
// #include <omp.h>


typedef double		 					 FT;
typedef COPT::Array<FT> 				Array;
typedef COPT::VectorBase<FT>			Vector;
typedef COPT::MatrixBase<FT>	 		Matrix;


// typedef double 			FT;
// typedef COPT::Vector<FT>			Vector;
// typedef COPT::Matrix<FT>			Matrix;


/*
 *	CBlas Test
 *
 */
// int main(int argc,char* argv[])
// {
// 	double *x= new double[5];
// 	x[0] = 1.0;
// 	x[1] = 2.0;
// 	x[2] = 3.0;
// 	x[3] = 4.0;
// 	x[4] = 5.0;
// 	double *y = new double[5];
// 	y[0] = 0.0;
// 	y[1] = 0.0;
// 	y[2] = 0.0;
// 	y[3] = 0.0;
// 	y[4] = 0.0;
// 	for ( int i = 0 ; i < 5 ; ++ i ){
// 		std::cout<<y[i]<<" ";
// 	}
// 	std::cout<<std::endl;
// 	blas::copt_blas_copy(5,x,1,y,1);
// 	for ( int i = 0 ; i < 5 ; ++ i ){
// 		std::cout<<y[i]<<" ";
// 	}
// 	std::cout<<std::endl;
// }

/*
 *		Array test
 */
// int main(int argc,char* argv[])
// {
// 	Array arr1(5);
// 	Array arr2(5);
// 	std::cout<<arr2<<std::endl;
// 	arr1[0] = 1.0;
// 	arr1[1] = 2.0;
// 	arr1[2] = 3.0;
// 	arr1[3] = 4.0;
// 	arr1[4] = 5.0;
// 	arr1.swap(arr2);
// 	std::cout<<"array 1 "<<arr1<<std::endl;
// 	std::cout<<"array 2 "<<arr2<<std::endl;
// 	arr1.copy(arr2);
// 	std::cout<<"array 1 "<<arr1<<std::endl;
// }

/*
 *			Vector test
 *
 */
// int main(int argc,char* argv[])
// {
// 	Vector vec1(5);
// 	Vector vec2(5);
// 	std::cout<<vec1<<std::endl;
// 	vec1[0] = 1.0;
// 	vec1[1] = 2.0;
// 	vec1[2] = 3.0;
// 	vec1[3] = 4.0;
// 	vec1[4] = 5.0;
// 	vec1.swap(vec2);
// 	std::cout<<"Vector 1 "<<vec1<<std::endl;
// 	std::cout<<"Vector 2 "<<vec2<<std::endl;
// 	vec1.copy(vec2);
// 	std::cout<<"Vector 1 "<<vec1<<std::endl;
// 	std::cout<<vec1.dot(vec2)<<std::endl;
// 	std::cout<<2.0*vec1<<std::endl;
// 	Matrix mat(2,5);
// 	mat(0,0) = 1.0;
// 	std::cout<<mat*vec1<<std::endl;
// }
/*LeastSquare test

*/
int main(int argc,char* argv[])
{
	Matrix A(4,2);
	A(0,0) = 1; A(0,1) = -1;
	A(1,0) = -1;A(1,1) = 1;
	A(2,0) = 2; A(2,1) = -2;
	A(3,0) = -3;A(3,1) = 1;
	Vector b(4);
	b[0] = 1;
	b[1] = 2;
	b[2] = 3;
	b[3] = 4;
	Vector x(2);
	COPT::LeastMeanSquareMethod(A,b,0.01,x);
	std::cout<<x<<std::endl;
	COPT::LeastSquareMethod(A,b,x);
	std::cout<<x<<std::endl;
}
/*
 *			Matrix test
 */
// int main(int argc,char* argv[])
// {
// 	Matrix iden = Matrix::identity(5,5);
// 	Vector vec(5);
// 	std::cout<<iden.solve(vec)<<std::endl;  
// }
/*
 *			Gradient method test
 */
// int main(int argc,char* argv[])
// { 
// 	COPT::RosenbrockFunction<Vector> func;
// 	Vector vec(2);
// 	vec[0] = 1.0;
// 	vec[1] = 0.0;
// 	std::cout<<func(vec)<<std::endl;
// 	std::cout<<func.gradient(vec)<<std::endl;
// 	std::cout<<func.hessian(vec)<<std::endl; 
// 	Vector x(2);
// 	x[0] = 1.0; 
// 	double tol = 1e-7;   
// 	int iter_max = 100; 
// 	// COPT::steepestDescentUsingBackTracking(
// 	// 	func,
// 	// 	0.7,
// 	// 	1e-4,
// 	// 	x,
// 	// 	tol,
// 	// 	iter_max
// 	// 	);
// 	COPT::newtonMethod(
// 		func,
// 		x,
// 		tol,
// 		iter_max);
// 	std::cout<<x<<" error "<<tol<<" iters "<<iter_max<<std::endl;
// }

// int main(int argc,char* argv[])
// {
// 	// #pragma omp parallel
// 	// printf("Hello from thread %d, nthreads %d\n",omp_get_thread_num(),omp_get_num_thread());
// 	/*
// 	 * conjugate gradient method test
// 	 */
// 	// Vector b(3),x(3);
// 	// Matrix A(3,3);
// 	// A(0,0) = 1;
// 	// A(1,1) = 1;
// 	// A(2,2) = 1;

// 	// b(0) = 1;b(1) = 1;b(2) = 1;

// 	// int iters = 100;
// 	// FT tol = 1e-6;
// 	// std::cout<<x<<std::endl;
// 	// COPT::conjugateGradientWithoutPrecondition(A,b,x,iters,tol);
// 	// std::cout<<"iteration used: "<<iters<<std::endl;
// 	// std::cout<<"final error estimated "<<tol<<std::endl;
// 	// std::cout<<x<<std::endl;


// 	/*
// 	 * steepest descent method test
// 	 */
// 	// COPT::TestQuadFunction<Vector> qf;
// 	// COPT::TestQuadFunctionWithDiff<Vector> qfwd;
	
// 	// // Vector gradient = qf.gradient(x0);
// 	// // Vector direction = -gradient;

// 	// clock_t one, two;
// 	// one = clock();
// 	// // COPT::findStepLengthBackTracking(qf,x0,gradient,direction,0.5,0.7,alpha,iters);
// 	// for ( int i = 0 ; i < 10000 ; ++ i ){
// 	// 	Vector x0(2);
// 	// 	FT tol = 1e-6;
// 	// 	int iters = 100;
// 	// 	// FT alpha = 1.0;
// 	// 	COPT::steepestDescentUsingBackTracking(qf,0.5,0.7,x0,tol,iters);
// 	// }
// 	// two= clock();
// 	// std::cout << "Runtime： " << (double)(two-one) * 1000.0 / CLOCKS_PER_SEC << " ms!" << std::endl;
// 	// one = clock();
// 	// // COPT::findStepLengthBackTracking(qf,x0,gradient,direction,0.5,0.7,alpha,iters);
// 	// for ( int i = 0 ; i < 10000 ; ++ i ){
// 	// 	Vector x0(2);
// 	// 	FT tol = 1e-6;
// 	// 	int iters = 100;
// 	// 	// FT alpha = 1.0;
// 	// 	COPT::steepestDescentUsingBackTracking(qfwd,0.5,0.7,x0,tol,iters);
// 	// }
// 	// two= clock();
// 	// std::cout << "Runtime： " << (double)(two-one) * 1000.0 / CLOCKS_PER_SEC << " ms!" << std::endl;
// 	// std::cout<<"final error is "<<tol <<std::endl;;
// 	// std::cout<<"final iteration is "<<iters<<std::endl;
// 	// std::cout<<"result is "<<x0<<std::endl;
 
// 	// double* data = new double[10];
// 	// for ( int i = 0 ; i < 10 ; ++ i )
// 	// 	data[i] =  1.0;
// 	// std::cout<<data[0]<<' '<<"1"<<std::endl;
// 	// COPT::Vector<double> vec(10);
// 	// COPT::Vector<double> v(vec);
// 	// std::cout<<vec.size()<<std::endl;
// 	// // vec[0] = 1.0;
// 	// std::cout<<"1"<<std::endl;
// 	// vec[0] = 2.0;
// 	// v[1] = 1.0;
// 	// COPT::Vector<double> vv = vec+v;
// 	// std::cout<<vv[0]<<std::endl;
// 	// std::cout<<vv(1)<<std::endl;
// 	// const COPT::Vector<double> vvv(vv);
// 	// std::cout<<vvv(2)<<std::endl;
// 	// COPT::Matrix<double> mat(2,3);
// 	// mat(0,1) = 2;
// 	// std::cout<<mat(0,1)<<std::endl;
// 	// std::cout<<vec[0]<<' '<<v[-1]<<" Hello World!"<<std::endl;
// 	// COPT::QuadFunction<double> cf(1.0,-2.0); 
// 	// COPT::RootSolver<COPT::ScalarFunction<double> > rs(cf);
// 	// std::cout<<rs.solve();
// 	// typedef COPT::Vector<double> Vector;

// 	// COPT::Vector<double> vec(2);
// 	// vec[0] = 1.0;
// 	// vec[1] = 2.0;
// 	// COPT::VectorCosineFunction<Vector> cf(vec);
// 	// Vector v(2);
// 	// v[0] = 1.57;
// 	// v[1] = 0;
// 	// std::cout<<cf(v)<<std::endl;
// 	// std::cout<<cf.gradient(v)<<std::endl; 
// 	// std::cout<<rs.
// 	// std::cout<<cf(1.57)<<std::endl;
// 	// COPT::ScalarDifferential<COPT::CosineFunction<double> > diff(cf);	
// 	// std::cout<<cf.diff(1.57)<<' '<<diff.diff(1.57)<<' '<<fabs(cf.diff(1.57)-diff.diff(1.57))<<std::endl;
// 	// COPT::SAFE_DELETE(data);
// }