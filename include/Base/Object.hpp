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

	enum ObjectType{

	};

public:

	typedef copt_object 					ObjectCategory;
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

	/** return the information about the COPT object */
	const std::string& info() {return __str;}

	/*			a clear function 
	 *			some classes of COPT might contain some pointer stuff which
	 *			has to be released when it is necessary.
	 */
	virtual void clear();
};

void COPTObject::clear()
{
}

}// End of namespace COPT
#endif // OBJECT_HPP__