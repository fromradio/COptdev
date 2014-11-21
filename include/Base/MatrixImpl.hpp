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

template<class ScalarType>
void MatrixBase<ScalarType>::blockFromMatrix(const MatrixBase& mat,const std::vector<size_t>& rownums,const std::vector<size_t>& colnums)
{
	this->resize(rownums.size(),colnums.size());
	for ( int r = 0 ; r < rownums.size() ; ++ r )
	{
		if(rownums[r]>=mat.rows())
			throw COException("Index out of range in matrix blocking!");
		for ( int c = 0 ; c < colnums.size() ; ++ c )
		{
			if(colnums[c]>=mat.cols())
				throw COException("Index out of range in matrix blocking!");
			this->operator()(r,c) = mat(rownums[r],colnums[c]);
		}
	}
}

template<class ScalarType>
void MatrixBase<ScalarType>::columnBlockFromMatrix(const MatrixBase& mat,const std::vector<size_t>& colnums)
{
	this->resize(mat.rows(),colnums.size());
	
	for ( int r = 0 ; r < mat.rows() ; ++ r )
	{
		for ( int c = 0 ; c < colnums.size() ; ++ c )
		{
			if(colnums[c]>=mat.cols())
				throw COException("Index out of range in matrix blocking!");
			this->operator()(r,c) = mat(r,colnums[c]);
		}
	}
}

template<class ScalarType>
void MatrixBase<ScalarType>::rowBlockFromMatrix(const MatrixBase& mat,const std::vector<size_t>& rownums)
{
	this->resize(rownums.size(),mat.cols());
	for ( int r = 0 ; r < rownums.size() ; ++ r )
	{
		if (rownums[r]>=mat.rows())
			throw COException("Index out of range in matrix blocking!");
		for ( int c = 0 ; c < mat.cols() ; ++ c )
		{
			this->operator()(r,c)=mat(rownums[r],c);
		}
	}
}

template<class ScalarType>
void MatrixBase<ScalarType>::combineAlongRow(const MatrixBase& m1,const MatrixBase& m2)
{
	MatrixBase::stCombineAlongRow(m1,m2,*this);
}

template<class ScalarType>
void MatrixBase<ScalarType>::combineAlongColumn(const MatrixBase& m1,const MatrixBase& m2)
{
	MatrixBase::stCombineAlongColumn(m1,m2,*this);
}

template<class ScalarType>
void MatrixBase<ScalarType>::stCombineAlongRow(const MatrixBase& m1,const MatrixBase& m2,MatrixBase& m)
{
	if(m1.cols()!=m2.cols())
		throw COException("Please make sure the column number of two matrices is the same before combination!");
	m.resize(m1.rows()+m2.rows(),m1.cols());
	for ( size_t i = 0 ; i < m1.cols() ; ++ i ){
		for ( size_t j = 0 ; j < m1.rows() ; ++ j ){
			m(j,i) = m1(j,i);
		}
		size_t n = m1.rows();
		for ( size_t j = 0 ; j < m2.rows() ; ++ j ){
			m(j+n,i) = m2(j,i);
		}
	}
}

template<class ScalarType>
void MatrixBase<ScalarType>::stCombineAlongColumn(const MatrixBase& m1,const MatrixBase& m2,MatrixBase& m)
{
	if(m1.rows()!=m2.rows())
		throw COException("Please make sure the row number of two matrices is the same before combination!");
	m.resize(m1.rows(),m1.cols()+m2.cols());
	for ( size_t i = 0 ; i < m1.rows() ; ++ i ){
		for (size_t j = 0 ; j < m1.cols() ; ++ j ){
			m(i,j) = m1(i,j);
		}
		size_t n = m1.cols();
		for (size_t j = 0 ; j < m2.cols() ; ++ j ){
			m(i,j+n)=m2(i,j);
		}
	}
}

}// End of namespace COPT

#endif