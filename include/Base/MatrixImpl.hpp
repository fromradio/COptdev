//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef MATRIX_IMPL_H
#define MATRIX_IMPL_H

namespace COPT
{


// template<class ScalarType>
// MatrixBase<ScalarType> operator* (const ScalarType s,const MatrixBase<ScalarType>& mat)
// {
	
// }
/*			create an identity matrix with all diagonal elements as s
 *			/param m: the number of rows
 *			/param n: the number of columns
 *			/param s: the value of scalar
 */
template<class ScalarType>
MatrixBase<ScalarType> MatrixBase<ScalarType>::identity(
	size_t m,
	size_t n,
	const ScalarType s)
{
	MatrixBase result(m,n);
	size_t min = std::min(m,n);
	for ( int i = 0 ; i < min ; ++ i )
		result(i,i) = s;
	return result;
}
}// End of namespace COPT

#endif