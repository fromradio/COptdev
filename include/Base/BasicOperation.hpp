#ifndef BASICOPERATION_H
#define BASICOPERATION_H

namespace COPT
{

/** COPTlong */
#ifndef COPTlong
#ifdef _WIN64
#define COPTlong __int64
#else
#define COPTlong long
#endif
#endif

#define longsize unsigned long

const double ZERO_THRESH = 1e-10;				// the threshold to judge whether a scalar is zero
const int    MAX_SEARCH = 10000;				// default maximum number of search

const double	DEFAULT_CONVERGE_ERROR = 1e-5; 	// default converge error
const double 	DEFAULT_STEP_FOR_DIFFERENTIAL = 1e-5;

const double INFTY = 1e10;

/*
 *				Judge that whether a scalar is zero
 */
template<class T>
inline bool IS_ZERO( T data )
{
	return fabs(data) < ZERO_THRESH ? true : false;
}

template<class T>
inline void SAFE_DELETE(T* value)
{
	if ( value ) { delete value;}
}

template<class T>
inline void SAFE_DELETE_ARRAY(T* array)
{
	if ( array ) { delete[] array;}
}

/*		Basic class for iterator
 */
class Iterator
{

};

/** 	base class without default constructor */
class nondefaultconstructor
{
	nondefaultconstructor();
};

/**		base class who is not copyable */
class noncopyable
{
	noncopyable(const noncopyable& );
};

} // End of namespace COPT
#endif