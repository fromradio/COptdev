#ifndef COEXCEPTION_H
#define COEXCEPTION_H

#include <exception>
#include <string>

// the base class of the exception that is used in the library
namespace COPT
{
class COException : public std::exception
{
public:
	COException(const char* str):__str("COpt Exception: "){ __str.append(str); }
	virtual const char* what() const throw()
	{
		return __str.c_str();
	}
	~COException() throw(){}
private:
	std::string			__str;
};
};

#endif