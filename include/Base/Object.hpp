// 		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef OBJECT_HPP__
#define OBJECT_HPP__

namespace COPT
{
class COPTObject
{
private:

	std::string 		__str;

public:

	COPTObject()
		:
		__str("This is an object of light-weight open source library COPT.")
	{}

	COPTObject(const std::string& str)
	{
		__str.append("This is an object of light-weight open source library COPT.");
		__str.append(str);
	}

	virtual ~COPTObject(){}
	const std::string& info() {return __str;}

};

}
#endif // OBJECT_HPP__