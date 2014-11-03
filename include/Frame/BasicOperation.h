#ifndef BASICOPERATION_H
#define BASICOPERATION_H

namespace COPT
{
template<class T>
inline void SAFE_DELETE(T* value)
{
	if ( value ) { delete value;}
}

template<class T>
inline void SAFE_DELETE_ARRAY(T* array)
{
	if ( array ) {delete[] array;}
}


const double ZERO_THRESH = 1e-10;
const int    MAX_SEARCH = 10000;
};
#endif