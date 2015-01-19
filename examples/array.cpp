#include "Core"

using namespace COPT;
typedef Array<double,int,Dynamic> DArray;

int main(int argc, char *argv[])
{
	/** a size-specified array */
	Array<double,int,3> a;
	a[0]=1.0;
	a[1]=2.0;
	a[2]=3.0;
	std::cout<<"a is "<<a<<std::endl;

	/** set the array */
	double s[3]{3.0,2.0,1.0};
	a.setArray(s);
	std::cout<<"a becomes "<<a<<std::endl;
}