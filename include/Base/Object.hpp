// 		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef OBJECT_HPP__
#define OBJECT_HPP__

namespace COPT
{
/*			The very base class of open source library 'COPT'. Actually
 *			there is nothing special in this class. Each class in COPT derives
 *			from 'COPTObject'. There are some virtual functions defined in 
 *			COPTObject. COPTObject stores two strings containing basic information of
 *			the class. One is the type of the class taking COPTObject as default.
 *			Another one is a small introduction of the class.
 *				virtual functions:
 *
 *				1. void clear();
 *					Clear the current object, nothing will be done if the derived
 *				class implement nonthing.
 *				2. const ostring& intro() const;
 *					return the brief introduction of the class.
 *	
 */
class COPTObject
{
private:

	/** a brief introduction of the object */
	ostring 		__intro;
	/** a string storing the type of the object */
	ostring			__type; 		

public:

	/** the category of the object */
	typedef copt_object 					ObjectCategory;

	/** default constructor of COPTObject */
	COPTObject(
		const ostring intro=ostring("This is an object of light-weight open source library COPT."),
		const ostring type=ostring("COPT object"))
		:
		__intro(intro),
		__type(type)
	{
	}

	virtual ~COPTObject(){}

	/** return the introduction about the COPT object */
	const ostring& intro() const;

	/** the type of the object */
	const ostring& type() const;

	/*			a clear function 
	 *			some classes of COPT might contain some pointer stuff which
	 *			has to be released when it is necessary.
	 */
	virtual void clear();
};

const ostring& COPTObject::intro() const
{
	return __intro;
}

const ostring& COPTObject::type() const
{
	return __type;
}

void COPTObject::clear()
{
}

}// End of namespace COPT
#endif // OBJECT_HPP__