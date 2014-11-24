//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef MATRIX_IMPL_H
#define MATRIX_IMPL_H

namespace COPT
{

template<class ScalarType>
MatrixBase<ScalarType>::MatrixBase()
	:
	Arr(),
	__rows(0),
	__cols(0)
{
}

template<class ScalarType>
MatrixBase<ScalarType>::MatrixBase(
	size_t m,
	size_t n,
	ScalarType* data)
	:
	Arr(m*n,data),
	__rows(m),
	__cols(n)
{
}

template<class ScalarType>
MatrixBase<ScalarType>::MatrixBase(
	const MatrixBase& mat)
	:
	Arr(mat.rows()*mat.cols(),mat.dataPtr()),
	__rows(mat.rows()),
	__cols(mat.cols())
{
}

template<class ScalarType>
MatrixBase<ScalarType>::~MatrixBase()
{
}

template<class ScalarType>
const size_t& MatrixBase<ScalarType>::rows() const
{
	return __rows;
}

template<class ScalarType>
const size_t& MatrixBase<ScalarType>::cols() const
{
	return __cols;
}

template<class ScalarType>
typename MatrixBase<ScalarType>::ScalarType& MatrixBase<ScalarType>::operator() (const int i,const int j)
{
	if(i<0||j<0)
		throw COException("MatrixBase error: index is less than zero!");
	else if (i>=__rows||j>=__cols)
		throw COException("MatrixBase error: index is out of range!");
	else{
		return this->operator[](j*__rows+i);
	}
}

template<class ScalarType>
const typename MatrixBase<ScalarType>::ScalarType& MatrixBase<ScalarType>::operator() ( const int i , const int j ) const
{
	return const_cast<MatrixBase&>(*this).operator()(i,j);
}

template<class ScalarType>
const typename MatrixBase<ScalarType>::ScalarType& MatrixBase<ScalarType>::data ( const int i ) const{
	if ( i < 0 )
		throw COException("MatrixBase error: index is less that zero!");
	else if ( i >= __rows*__cols )
		throw COException("MatrixBase error: index is out of range!");
	else
		return this->operator[](i);
}

template<class ScalarType>
VectorBase<ScalarType> MatrixBase<ScalarType>::row(const size_t num){
	if ( num >= __rows )
		throw COException("MatrixBase error: row index out of range!");
	else
		return VectorBase<ScalarType>(this->cols(),referred_array(),this->dataPtr()+num,this->rows());
}

template<class ScalarType>
VectorBase<ScalarType> MatrixBase<ScalarType>::col(const size_t num){
	if ( num >= __cols )
		throw COException("MatrixBase error: col index out of range!");
	else
		return VectorBase<ScalarType>(this->rows(),referred_array(),this->dataPtr()+num*this->rows(),1);
}

template<class ScalarType>
void MatrixBase<ScalarType>::resize(size_t m,size_t n)
{
	__rows = m;
	__cols = n;
	this->reset(m*n);
}

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

/**			Implementation of Triplet			*/

template<class ScalarType>
Triplet<ScalarType>::Triplet(
	const size_t r,
	const size_t c,
	const Scalar v)
	:
	__r(r),
	__c(c),
	__v(v)
{
}

template<class ScalarType>
Triplet<ScalarType>::~Triplet()
{
}

template<class ScalarType>
const size_t& Triplet<ScalarType>::rowIndex()
{
	return __r;
}

template<class ScalarType>
const size_t& Triplet<ScalarType>::columnIndex()
{
	return __c;
}

template<class ScalarType>
const ScalarType& Triplet<ScalarType>::value()
{
	return __v;
}

/**			Implementation of SpMatrixBase 		*/

template<class ScalarType>
SpMatrixBase<ScalarType>::SpMatrixBase()
	:
	__rows(0),
	__cols(0),
	__elesize(0),
	__rowind(NULL),
	__colptr(NULL),
	__vals(NULL)
{
}

template<class ScalarType>
SpMatrixBase<ScalarType>::SpMatrixBase(
	const size_t 					rows,
	const size_t 					cols,
	const size_t 					elesize,
	const size_t*					rowind,
	const size_t*					colptr,
	const ScalarType*				vals)
	:
	__rows(rows),
	__cols(cols),
	__elesize(elesize),
	__rowind(NULL),
	__colptr(NULL),
	__vals(NULL)
{
	__rowind = new size_t[__elesize];
	__vals = new ScalarType[__elesize];
	__colptr = new size_t[__cols+1];

	copt_blas_copy(__elesize,rowind,1,__rowind,1);
	copt_blas_copy(__elesize,vals,1,__vals,1);
	copt_blas_copy(__cols+1,colptr,1,__colptr,1);
}

template<class ScalarType>
SpMatrixBase<ScalarType>::~SpMatrixBase()
{
	if(__rowind)
		SAFE_DELETE_ARRAY(__rowind);
	if(__vals)
		SAFE_DELETE_ARRAY(__vals);
	if(__colptr)
		SAFE_DELETE_ARRAY(__colptr);
}

template<class ScalarType>
void SpMatrixBase<ScalarType>::setSparseMatrix(
	const size_t 					rows,
	const size_t 					cols,
	const size_t 					size,
	const size_t*					rowind,
	const size_t*			 		colptr,
	const ScalarType*			 	vals)
{
	clear();
	__rows = rows;
	__cols = cols;
	__elesize = size;
	
	__rows = new size_t[__elesize];
	__vals = new ScalarType[__elesize];
	__colptr = new size_t[__cols+1];

	copt_blas_copy(__elesize,rowind,1,__rowind,1);
	copt_blas_copy(__elesize,vals,1,__vals,1);
	copt_blas_copy(__cols+1,colptr,1,__colptr,1);
}

template<class ScalarType>
void SpMatrixBase<ScalarType>::setFromTriplets(
	const size_t rows,
	const size_t cols,
	const std::vector<Triplet>& triplets)
{
	clear();
	__rows = rows;
	__cols = cols;
}

template<class ScalarType>
void SpMatrixBase<ScalarType>::clear()
{
	__rows = 0;
	__cols = 0;
	__elesize = 0;

	if(__rowind)
	{
		SAFE_DELETE_ARRAY(__rowind);
		__rowind = NULL;
	}
	if(__vals)
	{
		SAFE_DELETE_ARRAY(__vals);
		__vals = NULL;
	}
	if(__colptr)
	{
		SAFE_DELETE_ARRAY(__colptr);
		__colptr = NULL;
	}
}

template<class ScalarType>
const ScalarType& SpMatrixBase<ScalarType>::operator()(
	const size_t i,
	const size_t j) const
{
	if(i>=__rows||j>=__cols)
		throw COException("Sparse Matrix error, index out of range!");
	size_t ip = __colptr[j],in=__colptr[j+1];
	for ( int i = ip ; i < in ; ++ i )
	{
		if( i == __rowind[i] )
			return __vals[i];
	}
	return __zero;
}

}// End of namespace COPT

#endif