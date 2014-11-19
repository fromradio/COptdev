//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef MATRIX_IMPL_H
#define MATRIX_IMPL_H

namespace COPT
{


template<class ScalarType>
void MatrixBase<ScalarType>::resize(size_t m,size_t n)
{
	__rows = m;
	__cols = n;
	this->reset(m*n);
}
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


template<class ScalarType>
void MatrixBase<ScalarType>::blockFromMatrix(const MatrixBase& mat,const std::set<size_t>& rownums,const std::set<size_t>& colnums)
{
	if (*rownums.rbegin()>=mat.rows()||*colnums.rbegin()>=mat.cols()){
		std::cerr<<*rownums.rbegin()<<' '<<*colnums.rbegin()<<std::endl;
		throw COException("Index out of range in matrix blocking!");
		return;
	}

	this->resize(rownums.size(),colnums.size());
	int r = 0,c = 0;
	for( std::set<size_t>::iterator riter = rownums.begin() ; riter != rownums.end() ; ++ riter)
	{
		c = 0;
		for ( std::set<size_t>::const_iterator citer = colnums.begin() ; citer != colnums.end() ; ++ citer )
		{
			this->operator()(r,c) = mat(*riter,*citer);
			++c;
		}
		++ r;
	}
}

template<class ScalarType>
void MatrixBase<ScalarType>::columnBlockFromMatrix(const MatrixBase& mat,const std::set<size_t>& colnums)
{
	if(*colnums.rbegin()>=mat.cols()){
		throw COException("Index out of range in matrix blocking!");
		return;
	}

	this->resize(mat.rows(),colnums.size());
	int c = 0;
	for ( std::set<size_t>::const_iterator citer = colnums.begin() ; citer != colnums.end() ; ++ citer )
	{
		for (int r = 0 ; r < mat.rows() ; ++ r )
		{
			this->operator()(r,c) = mat(r,*citer);
		}
		++ c;
	}
}

template<class ScalarType>
void MatrixBase<ScalarType>::rowBlockFromMatrix(const MatrixBase& mat,const std::set<size_t>& rownums)
{
	if(*rownums.rbegin()>=mat.rows()){
		throw COException("Index out of range in matrix blocking!");
		return;
	}

	this->resize(rownums.size(),mat.cols());
	int r = 0;
	for ( std::set<size_t>::const_iterator riter = rownums.begin() ; riter != rownums.end() ; ++ riter )
	{
		for ( int c = 0 ; c < mat.cols() ; ++ c )
		{
			this->operator()(r,c) = mat(*riter,c);
		}
		++ r;
	}
}
}// End of namespace COPT

#endif