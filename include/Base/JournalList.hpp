// Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
// Copyright (C) MathU

#ifndef JOURNAL_LIST_HPP__
#define JOURNAL_LIST_HPP__
namespace COPT
{
/*		Class JournalList is designed for recording
 *		the journal list of a solver containing the input,
 *		the iteration details and so on
 */
class JournalList
{
private:
	int 		__level;
public:
	JournalList();
};

/** 	The printer of COPT library
 *
 *
 */
class Printer
{
	int 		__print_level;

public:
	Printer()
		:
		__print_level(0)
	{
	}

	ostream& print(const COPTObject& obj , ostream& os );
};

}

#endif