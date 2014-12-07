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
	std::cout<<mat.transpose()*mat<<std::endl;
	std::cout<<mtm<<std::endl;

	clock_t t1,t2,at1,at2;
	at1 = clock();
	t1 = clock();
	Matrix m = Matrix::random(5000,4000);
	Vector v = Vector::random(4000);   	
	t2 = clock();
	std::cout<<"generation costs "<<(double)(t2-t1)/CLOCKS_PER_SEC<<std::endl;
	t1 = clock();
	Vector r = m*v;
	t2 = clock();
	std::cout<<" multiplicatioin costs "<<(double)(t2-t1)/CLOCKS_PER_SEC<<std::endl;
	// t1 = clock();
	m.mtm(mtm);
	t2 = clock();
	std::cout<<" level 3 operation costs "<<(double)(t2-t1)/CLOCKS_PER_SEC<<std::endl;
	m.transpose()*m;
	at2 = clock();
	std::cout<<" all costs "<<(double)(at2-at1)/CLOCKS_PER_SEC<<std::endl;
}