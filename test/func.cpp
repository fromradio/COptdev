/*
 * func.cpp
 *
 *  Created on: May 6, 2015
 *      Author: ruimin
 */

#include "Core"

int main(int argc, char *argv[])
{
//	COPT::RFunction<double,double> f = sin;
//	std::cout<<f(0.1)<<std::endl;
	std::function<double(double)> ff = sin;
	std::cout << ff(0.1) << std::endl;
	COPT::RFunction<double, double> f(sin);
	COPT::RFunction<double, double> fff([](double x)
	{	return x+0.1;});
	std::cout << ff(0.1) << std::endl;
	f = fff;
	std::cout << ff(0.1) << std::endl;
	std::list<double> l;
	l.resize(4);
	std::vector<double> input
	{ 1.0, 5.0, 2.0, 3.0 };
	for (auto iter = l.begin(); iter != l.end(); ++iter)
	{
		std::cout << *iter << std::endl;
	}
	fff.apply(input, l);
	for (auto iter = l.begin(); iter != l.end(); ++iter)
	{
		std::cout << *iter << std::endl;
	}
	std::cout<<fff.asApplyFunction()(l)<<std::endl;
}

