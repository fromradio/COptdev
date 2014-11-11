// 		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef VECTOR_IMPL_H
#define VECTOR_IMPL_H
/*
				The implementation of Vector
*/

namespace COPT
{
/*			Default construction
 *			a null vector with size zero is created
 */
template<class ScalarType>
VectorBase<ScalarType>::VectorBase()
	:Array<ScalarType>()
{
}
/*				
 *			multipy a vector's transpose
 *			/param vec: 		the vector which is transposed
 *			return value:		the result matrix
 */
template<class ScalarType>
MatrixBase<ScalarType> VectorBase<ScalarType>::mulTrans(const VectorBase<ScalarType>& vec) const
{
	size_t 					m = this->__size;
	size_t 					n = vec.size();
	MatrixBase<ScalarType> 	result(m,n);
	for ( int i = 0 ; i < m ; ++ i )
		for ( int j = 0 ; j < n ; ++ j )
			result(i,j) = this->operator[](i)*vec[j];
	return result;
}
}// End of namespace COPT


#endif