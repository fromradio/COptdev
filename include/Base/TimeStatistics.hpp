//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef TIME_STATISTICS_HPP__
#define TIME_STATISTICS_HPP__


namespace COPT
{

/** declaration */
class TimeStatistics;
class NoTimeStatistics;
class SolverTimeStatistics;

class TimeComputation
{
private:
	clock_t 		__begin;
	clock_t 		__end;
	bool 			__is_beginned;
public:
	TimeComputation();

	template<class Time>
	void timeBegin( const Time& );

	void timeBegin( const NoTimeStatistics& ,const no_time_stat_tag&);
	void timeBegin( const SolverTimeStatistics& , const solver_time_stat_tag& );

	template<class Time>
	void timeEnd( const Time& );
	void timeEnd( const NoTimeStatistics& ,const no_time_stat_tag& );
	void timeEnd( const SolverTimeStatistics& , const solver_time_stat_tag& );

	double currentTime();
};

TimeComputation::TimeComputation()
	:
	__is_beginned(false)
{
}

template<class Time>
void TimeComputation::timeBegin( const Time& stat )
{
	timeBegin(stat,typename Time::TimeCategory());
}

void TimeComputation::timeBegin( const NoTimeStatistics& , const no_time_stat_tag& )
{
}

void TimeComputation::timeBegin( const SolverTimeStatistics& , const solver_time_stat_tag& )
{
	if(__is_beginned)
		std::cerr<<"Time computation warning: last computation has not termined yet "<<std::endl;
	__begin = clock();
	__is_beginned = true;
}

template<class Time>
void TimeComputation::timeEnd( const Time& stat )
{
	timeEnd(stat,typename Time::TimeCategory());
}

void TimeComputation::timeEnd( const NoTimeStatistics& , const no_time_stat_tag & )
{
}

void TimeComputation::timeEnd(const SolverTimeStatistics& , const solver_time_stat_tag& )
{
	if(!__is_beginned)
		std::cerr<<"Time computation warning: time computation has not beginned yet "<<std::endl;
	__end = clock();
	__is_beginned = false;
}

double TimeComputation::currentTime()
{
	return (double)(__end-__begin)/CLOCKS_PER_SEC;
}


class TimeStatistics
{
private:
	TimeComputation 		__time_computation;
public:
	TimeComputation& timeComputer();

	virtual void printTimeInfo() = 0;
};

TimeComputation& TimeStatistics::timeComputer()
{
	return __time_computation;
}

class NoTimeStatistics
	:
	public TimeStatistics
{
private:
public:
	typedef 	no_time_stat_tag 		TimeCategory;

	void computationBegin();
	void computationEnd();

	void solvingBegin();
	void solvingEnd();

	void printTimeInfo();
};

void NoTimeStatistics::computationBegin()
{
	this->timeComputer().timeBegin(*this);
}

void NoTimeStatistics::computationEnd()
{
	this->timeComputer().timeEnd(*this);
}

void NoTimeStatistics::solvingBegin()
{
	this->timeComputer().timeBegin(*this);
}

void NoTimeStatistics::solvingEnd()
{
	this->timeComputer().timeEnd(*this);
}

void NoTimeStatistics::printTimeInfo()
{
}

class SolverTimeStatistics
	:
	public TimeStatistics
{
	/** the computation time */
	double 			__computation_time;
	/** the solve time */
	double 			__solving_time;
	/** whether it is set */
	bool 			__is_computed;

public:
	typedef 	solver_time_stat_tag 	TimeCategory;
	SolverTimeStatistics();

	void computationBegin();
	void computationEnd();

	void solvingBegin();
	void solvingEnd();

	void setComputationTime( const double t );
	void setSolvingTime(const double t);

	double computationTime() const;
	double solvingTime() const;

	void printTimeInfo();
};

SolverTimeStatistics::SolverTimeStatistics()
	:
	__computation_time(0.0),
	__solving_time(0.0),
	__is_computed(false)
{
}

void SolverTimeStatistics::computationBegin()
{
	this->timeComputer().timeBegin(*this);
}

void SolverTimeStatistics::computationEnd()
{
	this->timeComputer().timeEnd(*this);
	this->setComputationTime(this->timeComputer().currentTime());
}

void SolverTimeStatistics::solvingBegin()
{
	this->timeComputer().timeBegin(*this);
}

void SolverTimeStatistics::solvingEnd()
{
	this->timeComputer().timeEnd(*this);
	this->setSolvingTime(this->timeComputer().currentTime());
}

void SolverTimeStatistics::setComputationTime( const double t )
{
	__computation_time = t;
	__is_computed = true;
}

void SolverTimeStatistics::setSolvingTime( const double t )
{
	__solving_time = t;
}

void SolverTimeStatistics::printTimeInfo()
{
	std::cout<<"computation costs "<<__computation_time<<"s"<<std::endl;
	std::cout<<"solving costs "<<__solving_time<<"s"<<std::endl;
}

}

#endif