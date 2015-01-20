// This file is part of COPT, a light-weight C++ based optimization open source library
//
// Copyright (C) 2015 Ruimin Wang <ruimin.wang13@gmail.com>
// Copyright (C) 2015 MathU
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#ifndef MATRIX_IMPL_HPP__
#define MATRIX_IMPL_HPP__

namespace COPT
{

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::MatrixBase()
	:
	Array(),
	__rows(0),
	__cols(0),
	__sym(false),
	__trans(false)
{
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::MatrixBase(
	index m,
	index n,
	const scalar* data)
	:
	Array(m*n,data),
	__rows(m),
	__cols(n),
	__sym(false),
	__trans(false)
{
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::MatrixBase(
	const MatrixBase& mat)
	:
	Array(mat.rows()*mat.cols(),mat.dataPtr()),
	__rows(mat.rows()),
	__cols(mat.cols()),
	__sym(mat.isSymmetric()),
	__trans(mat.isTranspose())
{
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::~MatrixBase()
{
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
const index& MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::rows() const
{
	return __rows;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
bool MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::isRowDynamic() const
{
	return __is_row_dynamic;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
const index& MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::cols() const
{
	return __cols;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
bool MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::isColumnDynamic() const
{
	return __is_column_dynamic;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::scalar& MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::operator() (const index i,const index j)
{
	if(i<0||j<0)
		throw COException("MatrixBase error: index is less than zero!");
	else if (i>=__rows||j>=__cols)
		throw COException("MatrixBase error: index is out of range!");
	else if (__sym)
	{
		if ( i <= j )
			return this->operator[](j*__rows+i);
		else
			return this->operator[](i*__rows+j);
	}
	else if(__trans)
	{
		return this->operator[](i*__cols+j);
	}
	else{
		return this->operator[](j*__rows+i);
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
const typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::scalar& MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::operator() ( const index i , const index j ) const
{
	return const_cast<MatrixBase&>(*this).operator()(i,j);
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::set(const index i , const scalar value)
{
	if ( i < 0 )
		throw COException("MatrixBase error: index is less that zero!");
	else if ( i >= __rows*__cols )
		throw COException("MatrixBase error: index is out of range!");
	else
		this->operator[](i) = value;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
const typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::scalar& MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::data ( const index i ) const{
	if ( i < 0 )
		throw COException("MatrixBase error: index is less that zero!");
	else if ( i >= __rows*__cols )
		throw COException("MatrixBase error: index is out of range!");
	else
		return this->operator[](i);
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
VectorBase<scalar,index> MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::row(const index num){
	if ( num >= __rows || num < 0 )
		throw COException("MatrixBase error: row index out of range when getting a row of the matrix!");
	else if(__trans)
		return Vector(this->rows(),referred_array(),this->dataPtr()+num*this->rows(),1);
	else
		return Vector(this->cols(),referred_array(),this->dataPtr()+num,this->rows());
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
const VectorBase<scalar,index> MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::row(const index num )const {
	if ( num >= __rows || num < 0 )
		throw COException("MatrixBase error: row index out of range when getting a row of the matrix!");
	else if(__trans)
		return Vector(this->rows(),referred_array(),this->dataPtr()+num*this->rows(),1);
	else
		return Vector(this->cols(),referred_array(),this->dataPtr()+num,this->rows());
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
VectorBase<scalar,index> MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::col(const index num){
	if ( num >= __cols || num < 0 )
		throw COException("MatrixBase error: col index out of range when getting a column of the matrix!");
	else if(__trans)
		return Vector(this->cols(),referred_array(),this->dataPtr()+num,this->rows());
	else
		return Vector(this->rows(),referred_array(),this->dataPtr()+num*this->rows(),1);
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
const VectorBase<scalar,index> MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::col(const index num) const{
	if ( num >= __cols || num < 0 )
		throw COException("MatrixBase error: col index out of range when getting a column of the matrix!");
	else if(__trans)
		return Vector(this->cols(),referred_array(),this->dataPtr()+num,this->rows());
	else
		return Vector(this->rows(),referred_array(),this->dataPtr()+num*this->rows(),1);
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::setSymmetricFlag( bool sym )
{
	__sym = sym;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
bool MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::isSymmetric() const
{
	return __sym;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::setTransposeFlag( bool trans )
{
	__trans = trans;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
bool MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::isTranspose() const
{
	return __trans;
}



template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::resize(index m,index n)
{
	__rows = m;
	__cols = n;
	this->reset(m*n);
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>& MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::operator=(const MatrixBase& mat)
{
	if(__rows != mat.rows() || __cols != mat.cols() )
	{
		__rows = mat.rows();
		__cols = mat.cols();
		this->reset(__rows*__cols);
	}
	blas::copt_blas_copy(__rows*__cols,mat.dataPtr(),1,this->dataPtr(),1);
	__trans = mat.isTranspose();
	__sym = mat.isSymmetric();
	return *this;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime> MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::operator+(const MatrixBase& mat )
{
	if ( __rows != mat.rows() || __cols != mat.cols() )
	{
		throw COException("MatrixBase summation error: the size of two matrices are not consistent!");
	}
	MatrixBase result(*this);
	blas::copt_blas_axpy(this->size(),1.0,mat.dataPtr(),1,result.dataPtr(),1);
	return result;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime> MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::operator-(const MatrixBase& mat )
{
	if ( __rows != mat.rows() || __cols != mat.cols() )
	{
		throw COException("MatrixBase summation error: the size of two matrices are not consistent!");
	}
	MatrixBase result(*this);
	blas::copt_blas_axpy(this->size(),-1.0,mat.dataPtr(),1,result.dataPtr(),1);
	return result;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
VectorBase<scalar,index> MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::operator*(const VectorBase<scalar,index>& vec )const
{
	if ( __cols != vec.size() )
		throw COException("MatrixBase multiply error: the size of MatrixBase and vector are not consistent!");
	VectorBase<scalar,index> result(__rows);
	if (__sym)
		blas::copt_blas_symv(CblasColMajor,CblasUpper,__rows,1.0,this->dataPtr(),__rows,vec.dataPtr(),vec.interval(),0.0,result.dataPtr(),1);
	else if (__trans)
	{
		if( is_real<scalar>::value )
			blas::copt_blas_gemv(CblasColMajor,CblasTrans,__cols,__rows,1.0,this->dataPtr(),__cols,vec.dataPtr(),vec.interval(),0.0,result.dataPtr(),1);
		else if ( is_complex<scalar>::value )
			blas::copt_blas_gemv(CblasColMajor,CblasConjTrans,__cols,__rows,1.0,this->dataPtr(),__cols,vec.dataPtr(),vec.interval(),0.0,result.dataPtr(),1);
	}
	else
	{
		blas::copt_blas_gemv(CblasColMajor,CblasNoTrans,__rows,__cols,1.0,this->dataPtr(),__rows,vec.dataPtr(),vec.interval(),0.0,result.dataPtr(),1);
	}
	return result;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime> MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::operator*(const MatrixBase& mat )const
{
	if ( __cols != mat.rows() )
		throw COException("MatrixBase multiply error: the size of two matrices are not consistent!");
	MatrixBase result(__rows,mat.cols());
	if(__trans&&mat.isTranspose())
	{
		if( is_real<scalar>::value )
			blas::copt_blas_gemm(CblasColMajor,CblasTrans,CblasTrans,__rows,mat.cols(),__cols,1.0,this->dataPtr(),__cols,mat.dataPtr(),mat.cols(),0.0,result.dataPtr(),__rows);
		else if( is_complex<scalar>::value )
			blas::copt_blas_gemm(CblasColMajor,CblasConjTrans,CblasConjTrans,__rows,mat.cols(),__cols,1.0,this->dataPtr(),__cols,mat.dataPtr(),mat.cols(),0.0,result.dataPtr(),__rows);
	}
	else if( __trans )
	{
		std::cout<<"right"<<std::endl;
		if( is_real<scalar>::value )
			blas::copt_blas_gemm(CblasColMajor,CblasTrans,CblasNoTrans,__rows,mat.cols(),__cols,1.0,this->dataPtr(),__cols,mat.dataPtr(),__cols,0.0,result.dataPtr(),__rows);
		else if ( is_complex<scalar>::value ){
			std::cout<<"here"<<std::endl;
			blas::copt_blas_gemm(CblasColMajor,CblasConjTrans,CblasNoTrans,__rows,mat.cols(),__cols,1.0,this->dataPtr(),__cols,mat.dataPtr(),__cols,0.0,result.dataPtr(),__rows);
		}
	}
	else if (mat.isTranspose())
	{
		if( is_real<scalar>::value )
			blas::copt_blas_gemm(CblasColMajor,CblasNoTrans,CblasTrans,__rows,mat.cols(),__cols,1.0,this->dataPtr(),__rows,mat.dataPtr(),mat.cols(),0.0,result.dataPtr(),__rows);
		else if( is_complex<scalar>::value )
			blas::copt_blas_gemm(CblasColMajor,CblasNoTrans,CblasConjTrans,__rows,mat.cols(),__cols,1.0,this->dataPtr(),__rows,mat.dataPtr(),mat.cols(),0.0,result.dataPtr(),__rows);
	}
	else
		blas::copt_blas_gemm(CblasColMajor,CblasNoTrans,CblasNoTrans,__rows,mat.cols(),__cols,1.0,this->dataPtr(),__rows,mat.dataPtr(),mat.rows(),0.0,result.dataPtr(),__rows);
	return result;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime> MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::transpose() const
{
	MatrixBase result(__cols,__rows,this->dataPtr());
	result.setTransposeFlag(true);
	return result;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
VectorBase<scalar,index> MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::transMulti( const Vector& vec ) const
{
	if(__rows != vec.size() )
		throw COException("transpose multiplication error: the size of vector and matrix is not consistent!");
	Vector result(__cols);
	if(__trans)
		blas::copt_blas_gemv(CblasColMajor,CblasNoTrans,__cols,__rows,1.0,this->dataPtr(),__cols,vec.dataPtr(),vec.interval(),0.0,result.dataPtr(),1);
	else
		blas::copt_blas_gemv(CblasColMajor,CblasTrans,__rows,__cols,1.0,this->dataPtr(),__rows,vec.dataPtr(),vec.interval(),0.0,result.dataPtr(),1);
	return result;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime> MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::transMulti(const MatrixBase& mat ) const
{
	if(__rows != mat.rows() )
		throw COException("transpose multiplication error: the size of two matrices are not consistent!");
	MatrixBase result(__cols,mat.cols());
	if (__trans&&mat.isTranspose())
		blas::copt_blas_gemm(CblasColMajor,CblasNoTrans,CblasTrans,__cols,mat.cols(),__rows,1.0,this->dataPtr(),__cols,mat.dataPtr(),mat.cols(),0.0,result.dataPtr(),result.rows());
	else if(__trans)
		blas::copt_blas_gemm(CblasColMajor,CblasNoTrans,CblasNoTrans,__cols,mat.cols(),__rows,1.0,this->dataPtr(),__cols,mat.dataPtr(),mat.rows(),0.0,result.dataPtr(),result.rows());
	else if(mat.isTranspose())
		blas::copt_blas_gemm(CblasColMajor,CblasTrans,CblasTrans,__cols,mat.cols(),__rows,1.0,this->dataPtr(),__rows,mat.dataPtr(),mat.cols(),0.0,result.dataPtr(),result.rows());
	else
		blas::copt_blas_gemm(CblasColMajor,CblasTrans,CblasNoTrans,__cols,mat.cols(),__rows,1.0,
			this->dataPtr(),__rows,mat.dataPtr(),mat.rows(),0.0,result.dataPtr(),result.rows());
	return result;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime> MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::identity(
	index m,
	index n,
	const scalar s)
{
	MatrixBase result(m,n);
	index min = std::min(m,n);
	for ( index i = 0 ; i < min ; ++ i )
		result(i,i) = s;
	return result;
}


template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::blockFromMatrix(const MatrixBase& mat,const std::set<index>& rownums,const std::set<index>& colnums)
{
	if (*rownums.rbegin()>=mat.rows()||*colnums.rbegin()>=mat.cols()){
		std::cerr<<*rownums.rbegin()<<' '<<*colnums.rbegin()<<std::endl;
		throw COException("Index out of range in matrix blocking!");
		return;
	}

	this->resize(rownums.size(),colnums.size());
	index r = 0,c = 0;
	for( typename std::set<index>::iterator riter = rownums.begin() ; riter != rownums.end() ; ++ riter)
	{
		c = 0;
		for ( typename std::set<index>::const_iterator citer = colnums.begin() ; citer != colnums.end() ; ++ citer )
		{
			this->operator()(r,c) = mat(*riter,*citer);
			++c;
		}
		++ r;
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::columnBlockFromMatrix(const MatrixBase& mat,const std::set<index>& colnums)
{
	if(*colnums.rbegin()>=mat.cols()){
		throw COException("Index out of range in matrix blocking!");
		return;
	}

	this->resize(mat.rows(),colnums.size());
	index c = 0;
	for ( typename std::set<index>::const_iterator citer = colnums.begin() ; citer != colnums.end() ; ++ citer )
	{
		for (index r = 0 ; r < mat.rows() ; ++ r )
		{
			this->operator()(r,c) = mat(r,*citer);
		}
		++ c;
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::rowBlockFromMatrix(const MatrixBase& mat,const std::set<index>& rownums)
{
	if(*rownums.rbegin()>=mat.rows()){
		throw COException("Index out of range in matrix blocking!");
		return;
	}

	this->resize(rownums.size(),mat.cols());
	index r = 0;
	for ( typename std::set<index>::const_iterator riter = rownums.begin() ; riter != rownums.end() ; ++ riter )
	{
		for ( index c = 0 ; c < mat.cols() ; ++ c )
		{
			this->operator()(r,c) = mat(*riter,c);
		}
		++ r;
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::blockFromMatrix(const MatrixBase& mat,const std::vector<index>& rownums,const std::vector<index>& colnums)
{
	this->resize(rownums.size(),colnums.size());
	for ( index r = 0 ; r < rownums.size() ; ++ r )
	{
		if(rownums[r]>=mat.rows())
			throw COException("Index out of range in matrix blocking!");
		for ( index c = 0 ; c < colnums.size() ; ++ c )
		{
			if(colnums[c]>=mat.cols())
				throw COException("Index out of range in matrix blocking!");
			this->operator()(r,c) = mat(rownums[r],colnums[c]);
		}
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime> template<class InputIterator>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::blockFromMatrix(
	const MatrixBase& mat,
	const InputIterator& rowbegin,
	const InputIterator& rowend,
	const InputIterator& colbegin,
	const InputIterator& colend)
{
	// count the number of columns and rows at first
	index rownum = 0, colnum = 0;
	for (InputIterator iter = rowbegin ; iter != rowend ; ++ iter )
		++ rownum;
	for (InputIterator iter = colbegin ; iter != colend ; ++ iter )
		++ colnum;
	this->resize(rownum,colnum);
	index r = 0, c = 0;
	for ( InputIterator ri = rowbegin ; ri != rowend ; ++ ri )
	{
		c = 0;
		if ( *ri >= mat.rows () )
			throw COException("Index out of range in matrix blocking!");
		for ( InputIterator ci = colbegin ; ci != colend ; ++ ci )
		{
			if ( *ci >= mat.cols() )
				throw COException("Column index out of range in matrix blocking!");
			this->operator()(r,c) = mat(*ri,*ci);
			++ c;
		}
		++ r;
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::columnBlockFromMatrix(const MatrixBase& mat,const std::vector<index>& colnums)
{
	this->resize(mat.rows(),colnums.size());
	
	for ( index r = 0 ; r < mat.rows() ; ++ r )
	{
		for ( index c = 0 ; c < colnums.size() ; ++ c )
		{
			if(colnums[c]>=mat.cols())
				throw COException("Index out of range in matrix blocking!");
			this->operator()(r,c) = mat(r,colnums[c]);
		}
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime> template<class InputIterator>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::columnBlockFromMatrix(
	const MatrixBase& mat,
	const InputIterator& colbegin,
	const InputIterator& colend)
{
	// count the number of columns
	index colnum = 0;
	for ( InputIterator iter = colbegin ; iter != colend ; ++ iter )
		++ colnum;
	this->resize(mat.rows(),colnum);
	for ( index r = 0 ; r < mat.rows() ; ++ r )
	{
		index c = 0;
		for ( InputIterator ci = colbegin ; ci != colend ; ++ ci )
		{
			if( *ci >= mat.cols())
				throw COException("Column index out of range in matrix blocking!");
			this->operator()(r,c) = mat(r,*ci);
			++ c;
		}
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::rowBlockFromMatrix(const MatrixBase& mat,const std::vector<index>& rownums)
{
	this->resize(rownums.size(),mat.cols());
	for ( index r = 0 ; r < rownums.size() ; ++ r )
	{
		if (rownums[r]>=mat.rows())
			throw COException("Index out of range in matrix blocking!");
		for ( index c = 0 ; c < mat.cols() ; ++ c )
		{
			this->operator()(r,c)=mat(rownums[r],c);
		}
	}
}

template<class scalar, class index> template<class InputIterator>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::rowBlockFromMatrix(
	const MatrixBase& mat,
	const InputIterator& rowbegin,
	const InputIterator& rowend)
{
	// count the number of rows
	index rownum = 0;
	for ( InputIterator iter = rowbegin ; iter != rowend ; ++ iter )
		++ rownum;
	this->resize(rownum,mat.cols());
	index r = 0;
	for ( InputIterator ri = rowbegin ; ri != rowend ; ++ ri )
	{
		if ( *ri >= mat.rows() || *ri < 0 )
			throw COException("Index out of range in matrix blocking!");
		for ( index c = 0 ; c < mat.cols() ; ++ c )
		{
			this->operator()(r,c) = mat(*ri,c);
		}
		++ r;
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::combineAlongRow(const MatrixBase& m1,const MatrixBase& m2)
{
	MatrixBase::stCombineAlongRow(m1,m2,*this);
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::combineAlongColumn(const MatrixBase& m1,const MatrixBase& m2)
{
	MatrixBase::stCombineAlongColumn(m1,m2,*this);
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::stCombineAlongRow(const MatrixBase& m1,const MatrixBase& m2,MatrixBase& m)
{
	if(m1.cols()!=m2.cols())
		throw COException("Please make sure the column number of two matrices is the same before combination!");
	m.resize(m1.rows()+m2.rows(),m1.cols());
	for ( index i = 0 ; i < m1.cols() ; ++ i ){
		for ( index j = 0 ; j < m1.rows() ; ++ j ){
			m(j,i) = m1(j,i);
		}
		index n = m1.rows();
		for ( index j = 0 ; j < m2.rows() ; ++ j ){
			m(j+n,i) = m2(j,i);
		}
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::stCombineAlongColumn(const MatrixBase& m1,const MatrixBase& m2,MatrixBase& m)
{
	if(m1.rows()!=m2.rows())
		throw COException("Please make sure the row number of two matrices is the same before combination!");
	m.resize(m1.rows(),m1.cols()+m2.cols());
	for ( index i = 0 ; i < m1.rows() ; ++ i ){
		for (index j = 0 ; j < m1.cols() ; ++ j ){
			m(i,j) = m1(i,j);
		}
		index n = m1.cols();
		for (index j = 0 ; j < m2.cols() ; ++ j ){
			m(i,j+n)=m2(i,j);
		}
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::setRandom(const index rows,const index cols)
{
	if(rows<0||cols<0)
		throw COException("Please make sure that the number of row and column is bigger than zero!");
	std::uniform_real_distribution<typename get_pod_type<scalar>::type> unif(0.0,1.0);
	this->resize(rows,cols);
	for ( int i = 0 ; i < rows ; ++ i )
		for ( int j = 0 ; j < cols ; ++ j )
		{
			if(is_real<scalar>::value)
				ForceAssignment(unif(copt_rand_eng),this->operator()(i,j));
			else
				ForceAssignment(std::complex<typename get_pod_type<scalar>::type>(unif(copt_rand_eng),unif(copt_rand_eng)),this->operator()(i,j));
		}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime> MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::random(const index rows,const index cols)
{
	MatrixBase result;
	result.setRandom(rows,cols);
	return result;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::mtm( MatrixBase& mat ) const
{
	int m = this->rows();
	int n = this->cols();
	mat.resize(n,n);
	if ( is_real<scalar>::value )
		blas::copt_blas_syrk(CblasColMajor,CblasUpper,CblasTrans,n,m,1.0,this->dataPtr(),m,0.0,mat.dataPtr(),n);
	else
		blas::copt_blas_herk(CblasColMajor,CblasUpper,CblasConjTrans,n,m,1.0,this->dataPtr(),m,0.0,mat.dataPtr(),n);
	/** remember to assign the other half */
	for ( int i = 0 ; i < mat.rows() ; ++ i )
		for ( int j = 0 ; j < i ; ++ j )
			mat.operator[](i+j*mat.cols()) = mat.operator[](j+i*mat.cols());
	mat.setSymmetricFlag(true);
}

template<class Matrix>
class PartialEigenSolver;

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::podscalar MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::operationNorm() const
{
	MatrixBase mtm;
	this->mtm(mtm);
	PartialEigenSolver<MatrixBase> solver(mtm);
	podscalar e = solver.computeLargestEigenvalue();
	return std::sqrt(e);
}
///////////////End of implementation of 'MatrixBase'



}// End of namespace COPT

#endif