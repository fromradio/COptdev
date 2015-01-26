#include "Core"
#include <iostream>

int main(int argc, char* argv[])
{
	COPT::MSize<int> s{1,2};
	int k=2,l=4;
	std::cout<<s.m()<<' '<<s.n()<<std::endl;
	s = {k,l};
	std::cout<<s.m()<<' '<<s.n()<<std::endl;
}