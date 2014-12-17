#include <Header>

typedef COPT::KernelTrait<double,int> 		kernel;
typedef kernel::Vector 						Vector;
typedef COPT::LinearConstraint<kernel>		constraint;
int main(int argc,char *argv[])
{
	constraint con(2);
	Vector vec(2);
	vec(0)=1;
	std::cout<<con.feasible(vec)<<std::endl;
}