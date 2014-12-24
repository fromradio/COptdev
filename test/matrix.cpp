#include <Header>
 

typedef double		 					FT;
typedef COPT::KernelTrait<FT>				kernel;
typedef kernel::index 					Index;
// typedef COPT::Array<FT> 				Array;
typedef kernel::Vector			Vector;
typedef kernel::Matrix	 		Matrix;
// typedef COPT::KernelTrait<FT>			kernel;

int main(int argc,char* argv[])
{
	Matrix mat(1,2);
	mat(0,0) = 2;
	mat(0,1) = 1;
	Matrix mtm;
	mat.mtm(mtm);
	std::cout<<"m is "<<mat<<std::endl;
	std::cout<<mat.transpose()<<std::endl;   
	std::cout<<mat.transpose()*mat<<std::endl;
	std::cout<<mtm<<std::endl;
	std::cout<<mat.transMulti(mat)<<std::endl;

	Vector v(1);
	v(0) = 1.0;
	Matrix  m = mat.transpose();
	std::cout<<"test is "<<std::endl<<m.transMulti(mat.transpose());


	Vector vv(2);
	vv(0) = 1.0; vv(1) = 1.0;
	std::cout<<"squared norm is "<<vv.squaredNorm()<<std::endl;
	// std::cout<<m*vv <<std::endl;
	// std::cout<<mat*mtm<<std::endl;
	// std::cout<<mat.transpose()*mat.col(0)<<std::endl; 

	for ( int i = 0 ; i < 4 ; ++ i )
		std::cout<<mtm.dataPtr()[i]<<std::endl; 


	// clock_t t1,t2,at1,at2;
	// at1 = clock();
	// t1 = clock();
	// Matrix m = Matrix::random(5,2);
	// mat = m;
	// m = Matrix::random(5,2);
	// Vector v = Vector::random(4);
	// std::cout<<"m is "<<m<<std::endl;
	// std::cout<<"mat is "<<mat<<std::endl;
	// std::cout<<" m+mat is "<<m+mat<<std::endl;
	// std::cout<<" m-mat is "<<m-mat<<std::endl;
	// t2 = clock();
	// std::cout<<"generation costs "<<(double)(t2-t1)/CLOCKS_PER_SEC<<std::endl;
	// t1 = clock();
	// mat = m;
	// t2 = clock();
	// std::cout<<"assignment costs "<<(double)(t2-t1)/CLOCKS_PER_SEC<<std::endl;  
	// t1 = clock();
	// Vector r = m*v;
	// t2 = clock();
	// std::cout<<" multiplicatioin costs "<<(double)(t2-t1)/CLOCKS_PER_SEC<<std::endl;
	// // t1 = clock();
	// m.mtm(mtm);
	// t2 = clock();
	// std::cout<<" level 3 operation costs "<<(double)(t2-t1)/CLOCKS_PER_SEC<<std::endl;
	// m.transpose()*m;
	// at2 = clock();
	// std::cout<<" all costs "<<(double)(at2-at1)/CLOCKS_PER_SEC<<std::endl;
}