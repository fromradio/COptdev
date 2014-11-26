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
TripletBase<ScalarType>::TripletBase(
	const size_t r,
	const size_t c,
	const ScalarType v)
	:
	__r(r),
	__c(c),
	__v(v)
{
}

template<class ScalarType>
TripletBase<ScalarType>::~TripletBase()
{
}

template<class ScalarType>
const size_t& TripletBase<ScalarType>::rowIndex() const
{
	return __r;
}

template<class ScalarType>
const size_t& TripletBase<ScalarType>::columnIndex() const
{
	return __c;
}

template<class ScalarType>
const ScalarType& TripletBase<ScalarType>::value() const
{
	return __v;
}

/** 		Implementation of rowComparison 	*/
template<class Triplet>
bool rowComparison<Triplet>::operator() (const Triplet& t1,const Triplet& t2)
{
	return t1.rowIndex()<t2.rowIndex();
}

/**			Implementation of columnComparison	*/
template<class Triplet>
bool columnComparison<Triplet>::operator() (const Triplet& t1,const Triplet& t2)
{
	return t1.columnIndex()<t2.columnIndex();
}

/**			Implementation of SpMatrixBase 		*/

template<class ScalarType>
void SpMatrixBase<ScalarType>::judgeRationality()
{
	for ( int i = 0 ; i < __elesize ; ++ i )
	{
		if ( __rowind[i] >= __rows )
			throw COException("Sparse matrix not rational: row index out of range!");
	}
}

template<class ScalarType>
const typename SpMatrixBase<ScalarType>::ScalarType SpMatrixBase<ScalarType>::__zero = static_cast<ScalarType>(0.0);

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

	blas::copt_blas_copy(__elesize,rowind,1,__rowind,1);
	blas::copt_blas_copy(__elesize,vals,1,__vals,1);
	blas::copt_blas_copy(__cols+1,colptr,1,__colptr,1);

	judgeRationality();
}

template<class ScalarType>
SpMatrixBase<ScalarType>::SpMatrixBase(
	const SpMatrixBase& mat)
{
	setSparseMatrix(
		mat.rows(),
		mat.cols(),
		mat.elementSize(),
		mat.rowIndex(),
		mat.columnPointer(),
		mat.values());
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
	
	__rowind = new size_t[__elesize];
	__vals = new ScalarType[__elesize];
	__colptr = new size_t[__cols+1];

	blas::copt_blas_copy(__elesize,rowind,1,__rowind,1);
	blas::copt_blas_copy(__elesize,vals,1,__vals,1);
	blas::copt_blas_copy(__cols+1,colptr,1,__colptr,1);

	judgeRationality();
}

template<class ScalarType>
SpMatrixBase<ScalarType>& SpMatrixBase<ScalarType>::operator=(const SpMatrixBase& mat)
{
	setSparseMatrix(
		mat.rows(),
		mat.cols(),
		mat.elementSize(),
		mat.rowIndex(),
		mat.columnPointer(),
		mat.values());
	return *this;
}

template<class ScalarType>
void SpMatrixBase<ScalarType>::setFromTriplets(
	const size_t rows,
	const size_t cols,
	std::vector<Triplet>& triplets)
{
	clear();
	__rows = rows;
	__cols = cols;
	// sort the triplets according to the column index at first
	std::sort(triplets.begin(),triplets.end(),columnComparison<Triplet>());
	// compute how many elements there are in one column
	size_t colind = 0 , ip = 0;
	for ( size_t i = 0 ; i < triplets.size() ; ++ i )
	{
		if(triplets[i].columnIndex()>colind)
		{
			colind = triplets[i].columnIndex();
			std::sort(triplets.begin()+ip,triplets.begin()+i,rowComparison<Triplet>());
			ip = i;
		}
	}
	// last sort
	std::sort(triplets.begin()+ip,triplets.end(),rowComparison<Triplet>());
	std::vector<size_t> 	colcounts(__cols,0);
	std::list<size_t> 		rowinds;
	std::list<ScalarType>	vals;
	colind = 0 , ip = 0;
	size_t count = 0;
	for ( size_t i = 0 ; i < triplets.size() ; ++ i )
	{
		if(triplets[i].columnIndex()>colind)
		{
			colcounts[colind] = count;
			colind = triplets[i].columnIndex();
			count = 0;
		}
		
		if(i==0)
		{
			rowinds.push_back(triplets[i].rowIndex());
			vals.push_back(triplets[i].value());
			++count;
		}
		else if(triplets[i].rowIndex()==triplets[i-1].rowIndex())
		{
			vals.back() += triplets[i].value();
		}
		else
		{
			rowinds.push_back(triplets[i].rowIndex());
			vals.push_back(triplets[i].value());
			++count;
		}
	}
	// last column
	colcounts[colind] = count;

	__colptr = new size_t[__cols+1];
	size_t columncount = 0;
	for ( size_t i = 0 ; i < __cols ; ++ i )
	{
		__colptr[i] = columncount;
		columncount += colcounts[i];
	}
	__colptr[__cols] = columncount;

	__rowind = new size_t[rowinds.size()];
	size_t i = 0;
	for ( std::list<size_t>::iterator iter = rowinds.begin() ; iter != rowinds.end() ; ++ iter , ++ i )
		__rowind[i] = *iter;

	__vals = new ScalarType[vals.size()];
	i = 0;
	for ( typename std::list<ScalarType>::iterator iter = vals.begin() ; iter != vals.end() ; ++ iter , ++ i )
		__vals[i] = *iter;

	__elesize = rowinds.size();

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
const size_t& SpMatrixBase<ScalarType>::rows() const
{
	return __rows;
}

template<class ScalarType>
const size_t& SpMatrixBase<ScalarType>::cols() const
{
	return __cols;
}

template<class ScalarType>
const size_t& SpMatrixBase<ScalarType>::elementSize() const
{
	return __elesize;
}

template<class ScalarType>
const size_t* SpMatrixBase<ScalarType>::columnPointer() const
{
	return __colptr;
}

template<class ScalarType>
const size_t* SpMatrixBase<ScalarType>::rowIndex() const
{
	return __rowind;
}

template<class ScalarType>
const ScalarType* SpMatrixBase<ScalarType>::values() const
{
	return __vals;
}

template<class ScalarType>
void SpMatrixBase<ScalarType>::scale(const ScalarType s)
{
	for ( size_t i = 0 ; i < __elesize ; ++ i )
		__vals[i] *= s;
}

template<class ScalarType>
void SpMatrixBase<ScalarType>::neg()
{
	for ( size_t i = 0 ; i < __elesize ; ++ i )
		__vals[i] = -__vals[i];
}

template<class ScalarType>
const ScalarType& SpMatrixBase<ScalarType>::operator()(
	const size_t i,
	const size_t j) const
{
	if(i>=__rows||j>=__cols)
		throw COException("Sparse Matrix error, index out of range!");
	size_t ip = __colptr[j],in=__colptr[j+1];
	for ( size_t ind = ip ; ind < in ; ++ ind )
	{
		if( i == __rowind[ind] )
			return __vals[ind];
	}
	return __zero;
}

template<class ScalarType>
SpMatrixBase<ScalarType> SpMatrixBase<ScalarType>::operator*(const SpMatrixBase& mat) const
{
	if(__cols != mat.rows() )
		throw COException("Multiplication error: matrix size does not fit!");
	std::list<Triplet> tris;
	for ( size_t c = 0 ; c < __cols ; ++ c )
	{
		size_t ci = __colptr[c] , cn = __colptr[c+1];
		for ( size_t r = ci ; r < cn ; ++ r )
		{
			size_t rind = __rowind[r];
			// traverse mat
			for ( size_t mc = 0 ; mc < mat.cols() ; ++ mc )
			{
				size_t mci = mat.columnPointer()[mc], mcn = mat.columnPointer()[mc+1];
				for ( size_t mr = mci ; mr < mcn ; ++ mr )
				{
					if ( c == mat.rowIndex()[mr] )
						tris.push_back(Triplet(rind,mc,__vals[r]*mat.values()[mr]));
				}
			}
		}
	}
	std::vector<Triplet> triplets;
	triplets.reserve(tris.size());
	for ( typename std::list<Triplet>::iterator iter = tris.begin() ; iter != tris.end() ; ++ iter )
		triplets.push_back(*iter);
	SpMatrixBase result;
	result.setFromTriplets(__rows,mat.cols(),triplets);
	return result;
}


template<class ScalarType>
SpMatrixBase<ScalarType> SpMatrixBase<ScalarType>::operator+ ( const SpMatrixBase& mat ) const
{
	if( __cols != mat.cols() || __rows != mat.rows() )
		throw COException("Sparse matrix summation error: size does not fit!");
	std::list<size_t> inds;
	std::list<ScalarType> vals;
	size_t *colptr = new size_t[__cols+1];
	size_t count = 0;
	for ( size_t c = 0 ; c < __cols ; ++ c )
	{
		colptr[c] = count;
		size_t ci1 = __colptr[c],cn1 = __colptr[c+1],ci2 = mat.columnPointer()[c],cn2=mat.columnPointer()[c+1];
		size_t i1 = ci1,i2=ci2;
		while( i1<cn1 || i2<cn2 )
		{
			if( i1 < cn1 && i2 < cn2 )
			{
				if(__rowind[i1]==mat.rowIndex()[i2])
				{
					inds.push_back(__rowind[i1]);
					vals.push_back(__vals[i1]+mat.values()[i2]);
					++ i1;
					++ i2;
				}
				else if(__rowind[i1]<mat.rowIndex()[i2])
				{
					inds.push_back(__rowind[i1]);
					vals.push_back(__vals[i1]);
					++ i1;
				}
				else
				{
					inds.push_back(mat.rowIndex()[i2]);
					vals.push_back(mat.values()[i2]);
					++ i2;
				}
			}
			else if ( i1 < cn1 )
			{
				inds.push_back(__rowind[i1]);
				vals.push_back(__vals[i1]);
				++ i1;
			}
			else{
				inds.push_back(mat.rowIndex()[i2]);
				vals.push_back(mat.values()[i2]);
				++ i2;
			}
			++ count;
		}
	}
	colptr[__cols] = count;
	size_t *rowind = new size_t[inds.size()];
	size_t i = 0;
	for ( std::list<size_t>::iterator iter = inds.begin() ; iter != inds.end() ; ++ iter , ++i )
		rowind[i] = *iter;
	i = 0;
	ScalarType *vs = new ScalarType[vals.size()];
	for ( typename std::list<ScalarType>::iterator iter = vals.begin() ; iter != vals.end() ; ++ iter , ++ i )
		vs[i] = *iter;
	SpMatrixBase result;
	result.setSparseMatrix(__rows,__cols,inds.size(),rowind,colptr,vs);

	delete[]rowind;
	delete[]colptr;
	delete[]vs;
	return result;
}

template<class ScalarType>
SpMatrixBase<ScalarType> SpMatrixBase<ScalarType>::operator- ( const SpMatrixBase& mat ) const
{
	if( __cols != mat.cols() || __rows != mat.rows() )
		throw COException("Sparse matrix summation error: size does not fit!");
	std::list<size_t> inds;
	std::list<ScalarType> vals;
	size_t *colptr = new size_t[__cols+1];
	size_t count = 0;
	for ( size_t c = 0 ; c < __cols ; ++ c )
	{
		colptr[c] = count;
		size_t ci1 = __colptr[c],cn1 = __colptr[c+1],ci2 = mat.columnPointer()[c],cn2=mat.columnPointer()[c+1];
		size_t i1 = ci1,i2=ci2;
		while( i1<cn1 || i2<cn2 )
		{
			if( i1 < cn1 && i2 < cn2 )
			{
				if(__rowind[i1]==mat.rowIndex()[i2])
				{
					inds.push_back(__rowind[i1]);
					vals.push_back(__vals[i1]-mat.values()[i2]);
					++ i1;
					++ i2;
				}
				else if(__rowind[i1]<mat.rowIndex()[i2])
				{
					inds.push_back(__rowind[i1]);
					vals.push_back(__vals[i1]);
					++ i1;
				}
				else
				{
					inds.push_back(mat.rowIndex()[i2]);
					vals.push_back(-mat.values()[i2]);
					++ i2;
				}
			}
			else if ( i1 < cn1 )
			{
				inds.push_back(__rowind[i1]);
				vals.push_back(__vals[i1]);
				++ i1;
			}
			else{
				inds.push_back(mat.rowIndex()[i2]);
				vals.push_back(-mat.values()[i2]);
				++ i2;
			}
			++ count;
		}
	}
	colptr[__cols] = count;
	size_t *rowind = new size_t[inds.size()];
	size_t i = 0;
	for ( std::list<size_t>::iterator iter = inds.begin() ; iter != inds.end() ; ++ iter , ++i )
		rowind[i] = *iter;
	i = 0;
	ScalarType *vs = new ScalarType[vals.size()];
	for ( typename std::list<ScalarType>::iterator iter = vals.begin() ; iter != vals.end() ; ++ iter , ++ i )
		vs[i] = *iter;
	SpMatrixBase result;
	result.setSparseMatrix(__rows,__cols,inds.size(),rowind,colptr,vs);
	delete[]rowind;
	delete[]colptr;
	delete[]vs;
	return result;
}

template<class ScalarType>
SpMatrixBase<ScalarType> SpMatrixBase<ScalarType>::operator- () const
{
	SpMatrixBase result(*this);
	result.neg();
	return result;
}

template<class ScalarType>
typename SpMatrixBase<ScalarType>::Vector SpMatrixBase<ScalarType>::operator*(const Vector& vec) const
{
	if(__cols!=vec.size() )
		throw COException("Multiplication error: size does not fit!");
	Vector result(__rows);
	for ( size_t i = 0 ; i < __cols ; ++ i )
	{
		size_t ip = __colptr[i] , in = __colptr[i+1];
		for ( size_t r = ip ; r < in ; ++ r ){
			result[__rowind[r]] += __vals[r]*vec[i];
		}
	}
	return result;
}


template<class ScalarType,class T>
SpMatrixBase<ScalarType> operator* (const T s,const SpMatrixBase<ScalarType>& mat)
{
	return mat.operator*(static_cast<ScalarType>(s));
}

template<class ScalarType>
MatrixBase<ScalarType> SpMatrixBase<ScalarType>::toDenseMatrix() const
{
	MatrixBase<ScalarType> result(__rows,__cols);
	for ( size_t i = 0 ; i < __cols ; ++ i )
	{
		size_t ip = __colptr[i] , in = __colptr[i+1];
		for ( size_t r = ip ; r < in ; ++ r )
		{
			result(__rowind[r],i) = __vals[r];
		}
	}
	return result;
}

template<class ScalarType>
SpMatrixBase<ScalarType> SpMatrixBase<ScalarType>::operator* ( const ScalarType s ) const
{
	SpMatrixBase result(*this);
	result.scale(s);
	return result;
}

template<class ScalarType>
SpMatrixBase<ScalarType> operator*(const ScalarType s,const SpMatrixBase<ScalarType>& mat)
{
	return mat.operator*(s);
}

}// End of namespace COPT

#endif