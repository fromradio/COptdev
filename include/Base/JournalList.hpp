// This file is part of COPT, a light-weight C++ based optimization open source library
//
// Copyright (C) 2015 Ruimin Wang <ruimin.wang13@gmail.com>
// Copyright (C) 2015 MathU
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef JOURNAL_LIST_HPP__
#define JOURNAL_LIST_HPP__

namespace COPT
{
/*		Class JournalList is designed for recording
 *		the journal list of a solver containing the input,
 *		the iteration details and so on
 */
class JournalList
	:
	public COPTObject
{
protected:
	int 								__level;

	std::ostream* 						__os;
public:

	JournalList()
		:
		__os(NULL)
	{
		setPrintLevel(3);
	}
	virtual ~JournalList()
	{
		__os = NULL;
	}
	void setPrintLevel( const int level)
	{
		__level = level;
		if(__level==3)
			setStream(std::cout);
	};
	int printLevel() const
	{ 
		return __level;
	}
	void setStream(std::ostream&os)
	{
		__os = &os;
		__os->precision(8);
	}
	std::ostream& printStream()
	{
		return *__os;
	}
};



template<class Solver>
class SolverJournal
	:
	public JournalList,
	noncopyable
{
private:
	const Solver& 				__sol;

	SolverJournal();
public:
	SolverJournal( const Solver& sol );
	static inline std::ostream& solveBegin( const Solver& sol , std::ostream& , const int level );
	inline std::ostream& solveBegin();
	static inline std::ostream& solveEnd( const Solver& sol , std::ostream& , const int level  );
	inline std::ostream& solveEnd();
	static inline std::ostream& iterationBegin(const Solver& sol,std::ostream& , const int level  );
	inline std::ostream& iterationBegin();
	static inline std::ostream& iterationEnd(const Solver& sol,std::ostream& , const int level  );
	inline std::ostream& iterationEnd();
};

template<class Solver>
SolverJournal<Solver>::SolverJournal( const Solver& sol )
	:
	__sol(sol)
{
}

template<class Solver>
std::ostream& SolverJournal<Solver>::iterationBegin(const Solver& sol,std::ostream& os,const int level )
{
	return os;
}

template<class Solver>
std::ostream& SolverJournal<Solver>::iterationBegin( )
{
	return iterationBegin(this->__sol,*this->__os,this->__level);
}

template<class Solver>
std::ostream& SolverJournal<Solver>::iterationEnd(const Solver& sol,std::ostream& os,const int level )
{
	if(level==3)
		os<<sol.iterationNumber()<<"\t"<<sol.objective()<<"\t"<<sol.estimatedError()<<std::endl;
	return os;
}

template<class Solver>
std::ostream& SolverJournal<Solver>::iterationEnd()
{
	return iterationEnd(this->__sol,*this->__os,this->__level);
}

template<class Solver>
std::ostream& SolverJournal<Solver>::solveBegin( const Solver& sol , std::ostream& os,const int level )
{
	return os;
}

template<class Solver>
std::ostream& SolverJournal<Solver>::solveBegin()
{
	return solveBegin(this->__sol,*this->__os,this->__level);
}

template<class Solver>
std::ostream& SolverJournal<Solver>::solveEnd( const Solver& sol, std::ostream& os,const int level)
{
	if(level==3){
		os<<"final estimated error is "<<sol.estimatedError()<<std::endl;
		os<<"final result is "<<sol.result()<<std::endl;
	}
	return os;
}

template<class Solver>
std::ostream& SolverJournal<Solver>::solveEnd()
{
	return solveEnd(this->__sol,*this->__os,this->__level);
}

}

#endif