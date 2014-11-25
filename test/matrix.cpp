#include <Header>
 

typedef double		 					FT;
typedef COPT::Array<FT> 				Array;
typedef COPT::VectorBase<FT>			Vector;
typedef COPT::MatrixBase<FT>	 		Matrix;
typedef COPT::KernelTrait<FT>			kernel;

int main(int argc,char* argv[])
{
	Matrix mat(2,5);
	mat(1,0) = 2;
	mat(0,4) = 1;
	std::set<size_t> r;
	std::set<size_t> c;
	r.insert(0);
	// r.insert(1);
	c.insert(0);
	c.insert(4);
	Matrix m;
	m.blockFromMatrix(mat,r,c);
	std::cout<<"blocking test one:"<<std::endl<<m<<std::endl;
	m.rowBlockFromMatrix(mat,r);
	std::cout<<"blocking test two:"<<std::endl<<m<<std::endl;
	m.columnBlockFromMatrix(mat,c);
	std::cout<<"blocking test three :"<<std::endl<<m<<std::endl;

	Vector vec(5);
	vec[4] = 2.0;
	Vector v;
	v.blockFromVector(vec,c);
	std::cout<<vec.block(c)<<std::endl;
	std::cout<<v<<std::endl; 
}