#ifndef BASICOPERATION_H
#define BASICOPERATION_H

namespace COPT
{

#ifdef _WIN64
typedef __int64 		COPTlong;
#else
typedef long 			COPTlong;
#endif
typedef unsigned long 	longsize;
/** COPTlong */


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

/**		base class who is not copyable */
class noncopyable
{
	noncopyable(const noncopyable& );
	const noncopyable& operator=(const noncopyable&);
protected:
	noncopyable(){}
	~noncopyable(){}
};

} // End of namespace COPT
#endif