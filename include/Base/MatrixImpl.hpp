//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef MATRIX_IMPL_H
#define MATRIX_IMPL_H

namespace COPT
{

template<class ScalarType,class Size>
MatrixBase<ScalarType,Size>::MatrixBase()
	:
	Arr(),
	__rows(0),
	__cols(0)
{
}

template<class ScalarType,class Size>
MatrixBase<ScalarType,Size>::MatrixBase(
	Size m,
	Size n,
	ScalarType* data)
	:
	Arr(m*n,data),
	__rows(m),
	__cols(n)
{
}

template<class ScalarType,class Size>
MatrixBase<ScalarType,Size>::MatrixBase(
	const MatrixBase& mat)
	:
	Arr(mat.rows()*mat.cols(),mat.dataPtr()),
	__rows(mat.rows()),
	__cols(mat.cols())
{
}

template<class ScalarType,class Size>
MatrixBase<ScalarType,Size>::~MatrixBase()
{
}

template<class ScalarType,class Size>
const Size& MatrixBase<ScalarType,Size>::rows() const
{
	return __rows;
}

template<class ScalarType,class Size>
const Size& MatrixBase<ScalarType,Size>::cols() const
{
	return __cols;
}

template<class ScalarType,class Size>
typename MatrixBase<ScalarType,Size>::ScalarType& MatrixBase<ScalarType,Size>::operator() (const int i,const int j)
{
	if(i<0||j<0)
		throw COException("MatrixBase error: index is less than zero!");
	else if (i>=__rows||j>=__cols)
		throw COException("MatrixBase error: index is out of range!");
	else{
		return this->operator[](j*__rows+i);
	}
}

template<class ScalarType,class Size>
const typename MatrixBase<ScalarType,Size>::ScalarType& MatrixBase<ScalarType,Size>::operator() ( const int i , const int j ) const
{
	return const_cast<MatrixBase&>(*this).operator()(i,j);
}

template<class ScalarType,class Size>
const typename MatrixBase<ScalarType,Size>::ScalarType& MatrixBase<ScalarType,Size>::data ( const int i ) const{
	if ( i < 0 )
		throw COException("MatrixBase error: index is less that zero!");
	else if ( i >= __rows*__cols )
		throw COException("MatrixBase error: index is out of range!");
	else
		return this->operator[](i);
}

template<class ScalarType,class Size>
VectorBase<ScalarType,Size> MatrixBase<ScalarType,Size>::row(const Size num){
	if ( num >= __rows || num < 0 )
		throw COException("MatrixBase error: row index out of range!");
	else
		return Vector(this->cols(),referred_array(),this->dataPtr()+num,this->rows());
}

template<class ScalarType,class Size>
const VectorBase<ScalarType,Size> MatrixBase<ScalarType,Size>::row(const Size num )const {
	if ( num >= __rows || num < 0 )
		throw COException("MatrixBase error: row index out of range!");
	else
		return Vector(this->cols(),referred_array(),this->dataPtr()+num,this->rows());
}

template<class ScalarType,class Size>
VectorBase<ScalarType,Size> MatrixBase<ScalarType,Size>::col(const Size num){
	if ( num >= __cols || num < 0 )
		throw COException("MatrixBase error: col index out of range!");
	else
		return Vector(this->rows(),referred_array(),this->dataPtr()+num*this->rows(),1);
}

template<class ScalarType,class Size>
const VectorBase<ScalarType,Size> MatrixBase<ScalarType,Size>::col(const Size num) const{
	if ( num >= __cols || num < 0 )
		throw COException("MatrixBase error: col index out of range!");
	else
		return Vector(this->rows(),referred_array(),this->dataPtr()+num*this->rows(),1);
}

template<class ScalarType,class Size>
void MatrixBase<ScalarType,Size>::resize(Size m,Size n)
{
	__rows = m;
	__cols = n;
	this->reset(m*n);
}

template<class ScalarType,class Size>
MatrixBase<ScalarType,Size> MatrixBase<ScalarType,Size>::identity(
	Size m,
	Size n,
	const ScalarType s)
{
	MatrixBase result(m,n);
	Size min = std::min(m,n);
	for ( int i = 0 ; i < min ; ++ i )
		result(i,i) = s;
	return result;
}


template<class ScalarType,class Size>
void MatrixBase<ScalarType,Size>::blockFromMatrix(const MatrixBase& mat,const std::set<Size>& rownums,const std::set<Size>& colnums)
{
	if (*rownums.rbegin()>=mat.rows()||*colnums.rbegin()>=mat.cols()){
		std::cerr<<*rownums.rbegin()<<' '<<*colnums.rbegin()<<std::endl;
		throw COException("Index out of range in matrix blocking!");
		return;
	}

	this->resize(rownums.size(),colnums.size());
	int r = 0,c = 0;
	for( typename std::set<Size>::iterator riter = rownums.begin() ; riter != rownums.end() ; ++ riter)
	{
		c = 0;
		for ( typename std::set<Size>::const_iterator citer = colnums.begin() ; citer != colnums.end() ; ++ citer )
		{
			this->operator()(r,c) = mat(*riter,*citer);
			++c;
		}
		++ r;
	}
}

template<class ScalarType,class Size>
void MatrixBase<ScalarType,Size>::columnBlockFromMatrix(const MatrixBase& mat,const std::set<Size>& colnums)
{
	if(*colnums.rbegin()>=mat.cols()){
		throw COException("Index out of range in matrix blocking!");
		return;
	}

	this->resize(mat.rows(),colnums.size());
	int c = 0;
	for ( typename std::set<Size>::const_iterator citer = colnums.begin() ; citer != colnums.end() ; ++ citer )
	{
		for (int r = 0 ; r < mat.rows() ; ++ r )
		{
			this->operator()(r,c) = mat(r,*citer);
		}
		++ c;
	}
}

template<class ScalarType,class Size>
void MatrixBase<ScalarType,Size>::rowBlockFromMatrix(const MatrixBase& mat,const std::set<Size>& rownums)
{
	if(*rownums.rbegin()>=mat.rows()){
		throw COException("Index out of range in matrix blocking!");
		return;
	}

	this->resize(rownums.size(),mat.cols());
	int r = 0;
	for ( typename std::set<Size>::const_iterator riter = rownums.begin() ; riter != rownums.end() ; ++ riter )
	{
		for ( int c = 0 ; c < mat.cols() ; ++ c )
		{
			this->operator()(r,c) = mat(*riter,c);
		}
		++ r;
	}
}

template<class ScalarType,class Size>
void MatrixBase<ScalarType,Size>::blockFromMatrix(const MatrixBase& mat,const std::vector<Size>& rownums,const std::vector<Size>& colnums)
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

template<class ScalarType,class Size> template<class InputIterator>
void MatrixBase<ScalarType,Size>::blockFromMatrix(
	const MatrixBase& mat,
	const InputIterator& rowbegin,
	const InputIterator& rowend,
	const InputIterator& colbegin,
	const InputIterator& colend)
{
	// count the number of columns and rows at first
	Size rownum = 0, colnum = 0;
	for (InputIterator iter = rowbegin ; iter != rowend ; ++ iter )
		++ rownum;
	for (InputIterator iter = colbegin ; iter != colend ; ++ iter )
		++ colnum;
	this->resize(rownum,colnum);
	Size r = 0, c = 0;
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

template<class ScalarType,class Size>
void MatrixBase<ScalarType,Size>::columnBlockFromMatrix(const MatrixBase& mat,const std::vector<Size>& colnums)
{
	this->resize(mat.rows(),colnums.size());
	
	for ( Size r = 0 ; r < mat.rows() ; ++ r )
	{
		for ( Size c = 0 ; c < colnums.size() ; ++ c )
		{
			if(colnums[c]>=mat.cols())
				throw COException("Index out of range in matrix blocking!");
			this->operator()(r,c) = mat(r,colnums[c]);
		}
	}
}

template<class ScalarType,class Size> template<class InputIterator>
void MatrixBase<ScalarType,Size>::columnBlockFromMatrix(
	const MatrixBase& mat,
	const InputIterator& colbegin,
	const InputIterator& colend)
{
	// count the number of columns
	Size colnum = 0;
	for ( InputIterator iter = colbegin ; iter != colend ; ++ iter )
		++ colnum;
	this->resize(mat.rows(),colnum);
	for ( Size r = 0 ; r < mat.rows() ; ++ r )
	{
		Size c = 0;
		for ( InputIterator ci = colbegin ; ci != colend ; ++ ci )
		{
			if( *ci >= mat.cols())
				throw COException("Column index out of range in matrix blocking!");
			this->operator()(r,c) = mat(r,*ci);
			++ c;
		}
	}
}

template<class ScalarType,class Size>
void MatrixBase<ScalarType,Size>::rowBlockFromMatrix(const MatrixBase& mat,const std::vector<Size>& rownums)
{
	this->resize(rownums.size(),mat.cols());
	for ( Size r = 0 ; r < rownums.size() ; ++ r )
	{
		if (rownums[r]>=mat.rows())
			throw COException("Index out of range in matrix blocking!");
		for ( Size c = 0 ; c < mat.cols() ; ++ c )
		{
			this->operator()(r,c)=mat(rownums[r],c);
		}
	}
}

template<class ScalarType, class Size> template<class InputIterator>
void MatrixBase<ScalarType,Size>::rowBlockFromMatrix(
	const MatrixBase& mat,
	const InputIterator& rowbegin,
	const InputIterator& rowend)
{
	// count the number of rows
	Size rownum = 0;
	for ( InputIterator iter = rowbegin ; iter != rowend ; ++ iter )
		++ rownum;
	this->resize(rownum,mat.cols());
	Size r = 0;
	for ( InputIterator ri = rowbegin ; ri != rowend ; ++ ri )
	{
		if ( *ri >= mat.rows() || *ri < 0 )
			throw COException("Index out of range in matrix blocking!");
		for ( Size c = 0 ; c < mat.cols() ; ++ c )
		{
			this->operator()(r,c) = mat(*ri,c);
		}
	}
}

template<class ScalarType,class Size>
void MatrixBase<ScalarType,Size>::combineAlongRow(const MatrixBase& m1,const MatrixBase& m2)
{
	MatrixBase::stCombineAlongRow(m1,m2,*this);
}

template<class ScalarType,class Size>
void MatrixBase<ScalarType,Size>::combineAlongColumn(const MatrixBase& m1,const MatrixBase& m2)
{
	MatrixBase::stCombineAlongColumn(m1,m2,*this);
}

template<class ScalarType,class Size>
void MatrixBase<ScalarType,Size>::stCombineAlongRow(const MatrixBase& m1,const MatrixBase& m2,MatrixBase& m)
{
	if(m1.cols()!=m2.cols())
		throw COException("Please make sure the column number of two matrices is the same before combination!");
	m.resize(m1.rows()+m2.rows(),m1.cols());
	for ( Size i = 0 ; i < m1.cols() ; ++ i ){
		for ( Size j = 0 ; j < m1.rows() ; ++ j ){
			m(j,i) = m1(j,i);
		}
		Size n = m1.rows();
		for ( Size j = 0 ; j < m2.rows() ; ++ j ){
			m(j+n,i) = m2(j,i);
		}
	}
}

template<class ScalarType,class Size>
void MatrixBase<ScalarType,Size>::stCombineAlongColumn(const MatrixBase& m1,const MatrixBase& m2,MatrixBase& m)
{
	if(m1.rows()!=m2.rows())
		throw COException("Please make sure the row number of two matrices is the same before combination!");
	m.resize(m1.rows(),m1.cols()+m2.cols());
	for ( Size i = 0 ; i < m1.rows() ; ++ i ){
		for (Size j = 0 ; j < m1.cols() ; ++ j ){
			m(i,j) = m1(i,j);
		}
		Size n = m1.cols();
		for (Size j = 0 ; j < m2.cols() ; ++ j ){
			m(i,j+n)=m2(i,j);
		}
	}
}

/**			Implementation of Triplet			*/

template<class ScalarType,class Size>
TripletBase<ScalarType,Size>::TripletBase(
	const Size r,
	const Size c,
	const ScalarType v)
	:
	__r(r),
	__c(c),
	__v(v)
{
}

template<class ScalarType,class Size>
TripletBase<ScalarType,Size>::~TripletBase()
{
}

template<class ScalarType,class Size>
const Size& TripletBase<ScalarType,Size>::rowIndex() const
{
	return __r;
}

template<class ScalarType,class Size>
const Size& TripletBase<ScalarType,Size>::columnIndex() const
{
	return __c;
}

template<class ScalarType,class Size>
const ScalarType& TripletBase<ScalarType,Size>::value() const
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

template<class ScalarType,class Size>
void SpMatrixBase<ScalarType,Size>::judgeRationality()
{
	for ( int i = 0 ; i < __elesize ; ++ i )
	{
		if ( __rowind[i] >= __rows )
			throw COException("Sparse matrix not rational: row index out of range!");
	}
}

template<class ScalarType,class Size>
const typename SpMatrixBase<ScalarType,Size>::ScalarType SpMatrixBase<ScalarType,Size>::__zero = static_cast<ScalarType>(0.0);

template<class ScalarType,class Size>
SpMatrixBase<ScalarType,Size>::SpMatrixBase()
	:
	__rows(0),
	__cols(0),
	__elesize(0),
	__colptr(NULL),
	__rowind(NULL),
	__vals(NULL)
{
}

int umfpack_symbolic(
	int rows,int cols,
	const COPTlong* colptr,const COPTlong* rowind,const double* vals,void **Symbolic,
	const double* Control,double* Info);

template<class ScalarType,class Size>
SpMatrixBase<ScalarType,Size>::SpMatrixBase(
	const Size 						rows,
	const Size 						cols,
	const Size 						elesize,
	const Size*						colptr,
	const Size*						rowind,
	const ScalarType*				vals)
	:
	__rows(rows),
	__cols(cols),
	__elesize(elesize),
	__colptr(new Size[cols+1]),
	__rowind(new Size[elesize]),
	__vals(new ScalarType[elesize])
{

	blas::copt_blas_copy(__cols+1,colptr,1,__colptr,1);
	blas::copt_blas_copy(__elesize,rowind,1,__rowind,1);
	blas::copt_blas_copy(__elesize,vals,1,__vals,1);

	judgeRationality();
}

template<class ScalarType,class Size>
SpMatrixBase<ScalarType,Size>::SpMatrixBase(
	const SpMatrixBase& mat)
	:
	__colptr(NULL),
	__rowind(NULL),
	__vals(NULL)
{
	setSparseMatrix(
		mat.rows(),
		mat.cols(),
		mat.elementSize(),
		mat.columnPointer(),
		mat.rowIndex(),
		mat.values());
}

template<class ScalarType,class Size>
SpMatrixBase<ScalarType,Size>::~SpMatrixBase()
{
	if(__rowind)
		SAFE_DELETE_ARRAY(__rowind);
	if(__vals)
		SAFE_DELETE_ARRAY(__vals);
	if(__colptr)
		SAFE_DELETE_ARRAY(__colptr);
}

template<class ScalarType,class Size>
void SpMatrixBase<ScalarType,Size>::setSparseMatrix(
	const Size 					rows,
	const Size 					cols,
	const Size 					size,
	const Size*			 		colptr,
	const Size*					rowind,
	const ScalarType*			 	vals)
{
	clear();
	__rows = rows;
	__cols = cols;
	__elesize = size;
	
	__rowind = new Size[__elesize];
	__vals = new ScalarType[__elesize];
	__colptr = new Size[__cols+1];

	blas::copt_blas_copy(__elesize,rowind,1,__rowind,1);
	blas::copt_blas_copy(__elesize,vals,1,__vals,1);
	blas::copt_blas_copy(__cols+1,colptr,1,__colptr,1);

	judgeRationality();
}

template<class ScalarType,class Size>
SpMatrixBase<ScalarType,Size>& SpMatrixBase<ScalarType,Size>::operator=(const SpMatrixBase& mat)
{
	setSparseMatrix(
		mat.rows(),
		mat.cols(),
		mat.elementSize(),
		mat.columnPointer(),
		mat.rowIndex(),
		mat.values());
	return *this;
}

template<class ScalarType,class Size>
void SpMatrixBase<ScalarType,Size>::setFromTriplets(
	const Size rows,
	const Size cols,
	std::vector<Triplet>& triplets)
{
	clear();
	__rows = rows;
	__cols = cols;
	// sort the triplets according to the column index at first
	std::sort(triplets.begin(),triplets.end(),columnComparison<Triplet>());
	// compute how many elements there are in one column
	Size colind = 0 , ip = 0;
	for ( Size i = 0 ; i < triplets.size() ; ++ i )
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
	std::vector<Size> 	colcounts(__cols,0);
	std::list<Size> 		rowinds;
	std::list<ScalarType>	vals;
	colind = 0 , ip = 0;
	Size count = 0;
	for ( Size i = 0 ; i < triplets.size() ; ++ i )
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

	__colptr = new Size[__cols+1];
	Size columncount = 0;
	for ( Size i = 0 ; i < __cols ; ++ i )
	{
		__colptr[i] = columncount;
		columncount += colcounts[i];
	}
	__colptr[__cols] = columncount;

	__rowind = new Size[rowinds.size()];
	Size i = 0;
	for ( typename std::list<Size>::iterator iter = rowinds.begin() ; iter != rowinds.end() ; ++ iter , ++ i )
		__rowind[i] = *iter;

	__vals = new ScalarType[vals.size()];
	i = 0;
	for ( typename std::list<ScalarType>::iterator iter = vals.begin() ; iter != vals.end() ; ++ iter , ++ i )
		__vals[i] = *iter;

	__elesize = rowinds.size();
}

template<class ScalarType,class Size> template<class InputIterator>
void SpMatrixBase<ScalarType,Size>::setFromTriplets(
	const Size rows,
	const Size cols,
	const InputIterator& begin,
	const InputIterator& end)
{
	clear();
	__rows = rows;
	__cols = cols;
	// sort the triplets according to the column index at first
	std::sort(begin,end,columnComparison<Triplet>());
	// compute how many elements there are in one column
	Size colind = 0 , ip = 0;
	InputIterator previter = begin;
	for ( InputIterator iter = begin ; iter != end ; ++ iter )
	{
		if(iter->columnIndex()>colind)
		{
			colind = iter->columnIndex();
			std::sort(previter,iter,rowComparison<Triplet>());
			previter = iter;
		}
	}
	// last sort
	std::sort(previter,end,rowComparison<Triplet>());
	std::vector<Size> 		colcounts(__cols,0);
	std::list<Size> 		rowinds;
	std::list<ScalarType>	vals;
	colind = 0 , ip = 0;
	previter = begin;
	Size count = 0;
	for ( InputIterator iter = begin ; iter != end ; ++ iter )
	{
		if(iter->columnIndex()>colind)
		{
			colcounts[colind] = count;
			colind = iter->columnIndex();
			count = 0;
		}
		if(iter == begin)
		{
			rowinds.push_back(iter->rowIndex());
			vals.push_back(iter->value());
			++count;
		}
		else if(iter->rowIndex()==previter->rowIndex())
		{
			vals.back() += iter->value();
		}
		else
		{
			rowinds.push_back(iter->rowIndex());
			vals.push_back(iter->value());
			++count;
		}
		previter = iter;
	}
	// last column
	colcounts[colind] = count;

	__colptr = new Size[__cols+1];
	Size columncount = 0;
	for ( Size i = 0 ; i < __cols ; ++ i )
	{
		__colptr[i] = columncount;
		columncount += colcounts[i];
	}
	__colptr[__cols] = columncount;

	__rowind = new Size[rowinds.size()];
	Size i = 0;
	for ( typename std::list<Size>::iterator iter = rowinds.begin() ; iter != rowinds.end() ; ++ iter , ++ i )
		__rowind[i] = *iter;

	__vals = new ScalarType[vals.size()];
	i = 0;
	for ( typename std::list<ScalarType>::iterator iter = vals.begin() ; iter != vals.end() ; ++ iter , ++ i )
		__vals[i] = *iter;

	__elesize = rowinds.size();
}

template<class ScalarType,class Size>
void SpMatrixBase<ScalarType,Size>::clear()
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

template<class ScalarType,class Size>
const Size& SpMatrixBase<ScalarType,Size>::rows() const
{
	return __rows;
}

template<class ScalarType,class Size>
const Size& SpMatrixBase<ScalarType,Size>::cols() const
{
	return __cols;
}

template<class ScalarType,class Size>
const Size& SpMatrixBase<ScalarType,Size>::elementSize() const
{
	return __elesize;
}

template<class ScalarType,class Size>
const Size* SpMatrixBase<ScalarType,Size>::columnPointer() const
{
	return __colptr;
}

template<class ScalarType,class Size>
const Size* SpMatrixBase<ScalarType,Size>::rowIndex() const
{
	return __rowind;
}

template<class ScalarType,class Size>
const ScalarType* SpMatrixBase<ScalarType,Size>::values() const
{
	return __vals;
}

template<class ScalarType,class Size>
void SpMatrixBase<ScalarType,Size>::scale(const ScalarType s)
{
	for ( Size i = 0 ; i < __elesize ; ++ i )
		__vals[i] *= s;
}

template<class ScalarType,class Size>
void SpMatrixBase<ScalarType,Size>::neg()
{
	for ( Size i = 0 ; i < __elesize ; ++ i )
		__vals[i] = -__vals[i];
}

template<class ScalarType,class Size>
const ScalarType& SpMatrixBase<ScalarType,Size>::operator()(
	const Size i,
	const Size j) const
{
	if(i>=__rows||j>=__cols)
		throw COException("Sparse Matrix error, index out of range!");
	Size ip = __colptr[j],in=__colptr[j+1];
	for ( Size ind = ip ; ind < in ; ++ ind )
	{
		if( i == __rowind[ind] )
			return __vals[ind];
	}
	return __zero;
}

template<class ScalarType,class Size>
SpMatrixBase<ScalarType,Size> SpMatrixBase<ScalarType,Size>::operator*(const SpMatrixBase& mat) const
{
	if(__cols != mat.rows() )
		throw COException("Multiplication error: matrix size does not fit!");
	std::list<Triplet> tris;
	for ( Size c = 0 ; c < __cols ; ++ c )
	{
		Size ci = __colptr[c] , cn = __colptr[c+1];
		for ( Size r = ci ; r < cn ; ++ r )
		{
			Size rind = __rowind[r];
			// traverse mat
			for ( Size mc = 0 ; mc < mat.cols() ; ++ mc )
			{
				Size mci = mat.columnPointer()[mc], mcn = mat.columnPointer()[mc+1];
				for ( Size mr = mci ; mr < mcn ; ++ mr )
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


template<class ScalarType,class Size>
SpMatrixBase<ScalarType,Size> SpMatrixBase<ScalarType,Size>::operator+ ( const SpMatrixBase& mat ) const
{
	if( __cols != mat.cols() || __rows != mat.rows() )
		throw COException("Sparse matrix summation error: size does not fit!");
	std::list<Size> inds;
	std::list<ScalarType> vals;
	Size *colptr = new Size[__cols+1];
	Size count = 0;
	for ( Size c = 0 ; c < __cols ; ++ c )
	{
		colptr[c] = count;
		Size ci1 = __colptr[c],cn1 = __colptr[c+1],ci2 = mat.columnPointer()[c],cn2=mat.columnPointer()[c+1];
		Size i1 = ci1,i2=ci2;
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
	Size *rowind = new Size[inds.size()];
	Size i = 0;
	for ( typename std::list<Size>::iterator iter = inds.begin() ; iter != inds.end() ; ++ iter , ++i )
		rowind[i] = *iter;
	i = 0;
	ScalarType *vs = new ScalarType[vals.size()];
	for ( typename std::list<ScalarType>::iterator iter = vals.begin() ; iter != vals.end() ; ++ iter , ++ i )
		vs[i] = *iter;
	SpMatrixBase result;
	result.setSparseMatrix(__rows,__cols,inds.size(),colptr,rowind,vs);

	delete[]rowind;
	delete[]colptr;
	delete[]vs;
	return result;
}

template<class ScalarType,class Size>
SpMatrixBase<ScalarType,Size> SpMatrixBase<ScalarType,Size>::operator- ( const SpMatrixBase& mat ) const
{
	if( __cols != mat.cols() || __rows != mat.rows() )
		throw COException("Sparse matrix summation error: size does not fit!");
	std::list<Size> inds;
	std::list<ScalarType> vals;
	Size *colptr = new Size[__cols+1];
	Size count = 0;
	for ( Size c = 0 ; c < __cols ; ++ c )
	{
		colptr[c] = count;
		Size ci1 = __colptr[c],cn1 = __colptr[c+1],ci2 = mat.columnPointer()[c],cn2=mat.columnPointer()[c+1];
		Size i1 = ci1,i2=ci2;
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
	Size *rowind = new Size[inds.size()];
	Size i = 0;
	for ( typename std::list<Size>::iterator iter = inds.begin() ; iter != inds.end() ; ++ iter , ++i )
		rowind[i] = *iter;
	i = 0;
	ScalarType *vs = new ScalarType[vals.size()];
	for ( typename std::list<ScalarType>::iterator iter = vals.begin() ; iter != vals.end() ; ++ iter , ++ i )
		vs[i] = *iter;
	SpMatrixBase result;
	result.setSparseMatrix(__rows,__cols,inds.size(),colptr,rowind,vs);
	delete[]rowind;
	delete[]colptr;
	delete[]vs;
	return result;
}

template<class ScalarType,class Size>
SpMatrixBase<ScalarType,Size> SpMatrixBase<ScalarType,Size>::operator- () const
{
	SpMatrixBase result(*this);
	result.neg();
	return result;
}

template<class ScalarType,class Size>
typename SpMatrixBase<ScalarType,Size>::Vector SpMatrixBase<ScalarType,Size>::operator*(const Vector& vec) const
{
	if(__cols!=vec.size() )
		throw COException("Multiplication error: size does not fit!");
	Vector result(__rows);
	for ( Size i = 0 ; i < __cols ; ++ i )
	{
		Size ip = __colptr[i] , in = __colptr[i+1];
		for ( Size r = ip ; r < in ; ++ r ){
			result[__rowind[r]] += __vals[r]*vec[i];
		}
	}
	return result;
}


template<class ScalarType,class Size, class T>
SpMatrixBase<ScalarType,Size> operator* (const T s,const SpMatrixBase<ScalarType,Size>& mat)
{
	return mat.operator*(static_cast<ScalarType>(s));
}

template<class ScalarType,class Size>
MatrixBase<ScalarType,Size> SpMatrixBase<ScalarType,Size>::toDenseMatrix() const
{
	MatrixBase<ScalarType,Size> result(__rows,__cols);
	for ( Size i = 0 ; i < __cols ; ++ i )
	{
		Size ip = __colptr[i] , in = __colptr[i+1];
		for ( Size r = ip ; r < in ; ++ r )
		{
			result(__rowind[r],i) = __vals[r];
		}
	}
	return result;
}

template<class ScalarType,class Size>
VectorBase<ScalarType,Size> SpMatrixBase<ScalarType,Size>::solve(const VectorBase<ScalarType,Size>& vec)
{
	return UMFLinearSolver<SpMatrixBase>(*this).solve(vec);
}

template<class ScalarType,class Size>
SpMatrixBase<ScalarType,Size> SpMatrixBase<ScalarType,Size>::operator* ( const ScalarType s ) const
{
	SpMatrixBase result(*this);
	result.scale(s);
	return result;
}

template<class ScalarType,class Size>
SpMatrixBase<ScalarType,Size> operator*(const ScalarType s,const SpMatrixBase<ScalarType,Size>& mat)
{
	return mat.operator*(s);
}

}// End of namespace COPT

#endif