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

/** 	The printer of COPT library
 *
 *
 */
class Printer
	:
	public COPTObject
{
	int 		__print_level;

	/** data objects */
	static inline std::ostream& printObjectTag( const data_object& , std::ostream& os);
	static inline std::ostream& printObjectTag( const matrix_object& , std::ostream& os);
	static inline std::ostream& printObjectTag( const vector_object& , std::ostream& os );

	/** problems */
	static inline std::ostream& printObjectTag( const lasso_problem& , std::ostream& os );

	/** solvers */
	static inline std::ostream& printObjectTag( const proximal_solver& , std::ostream& os );
	static inline std::ostream& printObjectTag( const admm_solver& , std::ostream& os );
	static inline std::ostream& printObjectTag( const fista_solver& , std::ostream& os );

public:
	Printer()
		:
		__print_level(0)
	{
	}

	std::ostream& print(const COPTObject& obj , std::ostream& os );

	template<class Object>
	static inline std::ostream& printType( const Object& , std::ostream& os );
};

template<class Object>
std::ostream& Printer::printType( const Object& obj , std::ostream& os )
{
	return printObjectTag(typename Object::ObjectCategory(),os);
}

std::ostream& Printer::printObjectTag( const data_object& , std::ostream& os )
{
	os<<"-----------------This is a data object in COPT------------------"<<std::endl;
	return os;
}

std::ostream& Printer::printObjectTag( const matrix_object& , std::ostream& os )
{
	os<<"--------------------This is a matrix in COPT--------------------"<<std::endl;
	return os;
}

std::ostream& Printer::printObjectTag( const vector_object& , std::ostream& os )
{
	os<<"---------------------This is a vector in COPT-------------------"<<std::endl;
	return os;
}

std::ostream& Printer::printObjectTag( const lasso_problem& , std::ostream& os )
{
	os<<"------------------This is a lasso problem in COPT---------------"<<std::endl;
	return os;
}

std::ostream& Printer::printObjectTag( const proximal_solver& , std::ostream& os )
{
	os<<"-----------------This is a proximal problem in COPT-------------"<<std::endl;
	return os;
}

std::ostream& Printer::printObjectTag( const admm_solver& , std::ostream& os )
{
	os<<"------------------This is an ADMM solver in COPT----------------"<<std::endl;
	return os;
}

std::ostream& Printer::printObjectTag( const fista_solver& , std::ostream& os )
{
	os<<"------------------This is an FISTA solver in COPT---------------"<<std::endl;
	return os;
}

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
	if(level==3)
		os<<"final estimated error is "<<sol.estimatedError()<<std::endl;
	return os;
}

template<class Solver>
std::ostream& SolverJournal<Solver>::solveEnd()
{
	return solveEnd(this->__sol,*this->__os,this->__level);
}

}

#endif