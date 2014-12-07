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
	// std::set<Index> r;
	// std::set<Index> c;
	// r.insert(0);
	// // r.insert(1);
	// c.insert(0);
	// c.insert(4);
	// Matrix m;
	// m.blockFromMatrix(mat,r.begin(),r.end(),c.begin(),c.end());
	// std::cout<<"blocking test one:"<<std::endl<<m<<std::endl;
	// m.rowBlockFromMatrix(mat,r.begin(),r.end());
	// std::cout<<"blocking test two:"<<std::endl<<m<<std::endl;
	// m.columnBlockFromMatrix(mat,c.begin(),c.end());
	// std::cout<<"blocking test three :"<<std::endl<<m<<std::endl;

	// Vector vec(5);
	// vec[4] = 2.0;
	// Vector v;
	// v.blockFromVector(vec,c);
	// std::cout<<vec.block(c)<<std::endl;
    
	Matrix mtm;
	mat.mtm(mtm);
	std::cout<<"m is "<<mat<<std::endl;
	std::cout<<mat.transpose()*mat<<std::endl;
	std::cout<<mtm<<std::endl;      
}