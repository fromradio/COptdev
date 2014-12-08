//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef MATRIX_IMPL_H
#define MATRIX_IMPL_H

namespace COPT
{

template<class scalar,class index>
MatrixBase<scalar,index>::MatrixBase()
	:
	Arr(),
	__rows(0),
	__cols(0),
	__sym(false)
{
}

template<class scalar,class index>
MatrixBase<scalar,index>::MatrixBase(
	index m,
	index n,
	scalar* data)
	:
	Arr(m*n,data),
	__rows(m),
	__cols(n),
	__sym(false)
{
}

template<class scalar,class index>
MatrixBase<scalar,index>::MatrixBase(
	const MatrixBase& mat)
	:
	Arr(mat.rows()*mat.cols(),mat.dataPtr()),
	__rows(mat.rows()),
	__cols(mat.cols()),
	__sym(false)
{
}

template<class scalar,class index>
MatrixBase<scalar,index>::~MatrixBase()
{
}

template<class scalar,class index>
const index& MatrixBase<scalar,index>::rows() const
{
	return __rows;
}

template<class scalar,class index>
const index& MatrixBase<scalar,index>::cols() const
{
	return __cols;
}

template<class scalar,class index>
typename MatrixBase<scalar,index>::scalar& MatrixBase<scalar,index>::operator() (const index i,const index j)
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
	else{
		return this->operator[](j*__rows+i);
	}
}

template<class scalar,class index>
const typename MatrixBase<scalar,index>::scalar& MatrixBase<scalar,index>::operator() ( const index i , const index j ) const
{
	return const_cast<MatrixBase&>(*this).operator()(i,j);
}

template<class scalar,class index>
void MatrixBase<scalar,index>::set(const index i , const scalar value)
{
	if ( i < 0 )
		throw COException("MatrixBase error: index is less that zero!");
	else if ( i >= __rows*__cols )
		throw COException("MatrixBase error: index is out of range!");
	else
		this->operator[](i) = value;
}

template<class scalar,class index>
const typename MatrixBase<scalar,index>::scalar& MatrixBase<scalar,index>::data ( const index i ) const{
	if ( i < 0 )
		throw COException("MatrixBase error: index is less that zero!");
	else if ( i >= __rows*__cols )
		throw COException("MatrixBase error: index is out of range!");
	else
		return this->operator[](i);
}

template<class scalar,class index>
VectorBase<scalar,index> MatrixBase<scalar,index>::row(const index num){
	if ( num >= __rows || num < 0 )
		throw COException("MatrixBase error: row index out of range!");
	else
		return Vector(this->cols(),referred_array(),this->dataPtr()+num,this->rows());
}

template<class scalar,class index>
const VectorBase<scalar,index> MatrixBase<scalar,index>::row(const index num )const {
	if ( num >= __rows || num < 0 )
		throw COException("MatrixBase error: row index out of range!");
	else
		return Vector(this->cols(),referred_array(),this->dataPtr()+num,this->rows());
}

template<class scalar,class index>
VectorBase<scalar,index> MatrixBase<scalar,index>::col(const index num){
	if ( num >= __cols || num < 0 )
		throw COException("MatrixBase error: col index out of range!");
	else
		return Vector(this->rows(),referred_array(),this->dataPtr()+num*this->rows(),1);
}

template<class scalar,class index>
const VectorBase<scalar,index> MatrixBase<scalar,index>::col(const index num) const{
	if ( num >= __cols || num < 0 )
		throw COException("MatrixBase error: col index out of range!");
	else
		return Vector(this->rows(),referred_array(),this->dataPtr()+num*this->rows(),1);
}

template<class scalar,class index>
void MatrixBase<scalar,index>::setSymmetricFlag( bool sym )
{
	__sym = sym;
}

template<class scalar,class index>
bool MatrixBase<scalar,index>::isSymmetric() const
{
	return __sym;
}

template<class scalar,class index>
void MatrixBase<scalar,index>::resize(index m,index n)
{
	__rows = m;
	__cols = n;
	this->reset(m*n);
}

template<class scalar,class index>
MatrixBase<scalar,index>& MatrixBase<scalar,index>::operator=(const MatrixBase& mat)
{
	if(__rows != mat.rows() || __cols != mat.cols() )
	{
		__rows = mat.rows();
		__cols = mat.cols();
		this->reset(__rows*__cols);
	}
	blas::copt_blas_copy(__rows*__cols,mat.dataPtr(),1,this->dataPtr(),1);
	return *this;
}

template<class scalar,class index>
MatrixBase<scalar,index> MatrixBase<scalar,index>::operator+(const MatrixBase& mat )
{
	if ( __rows != mat.rows() || __cols != mat.cols() )
	{
		throw COException("MatrixBase summation error: the size of two matrices are not consistent!");
	}
	MatrixBase result(*this);
	blas::copt_blas_axpy(this->size(),1.0,mat.dataPtr(),1,result.dataPtr(),1);
	return result;
}

template<class scalar,class index>
MatrixBase<scalar,index> MatrixBase<scalar,index>::operator-(const MatrixBase& mat )
{
	if ( __rows != mat.rows() || __cols != mat.cols() )
	{
		throw COException("MatrixBase summation error: the size of two matrices are not consistent!");
	}
	MatrixBase result(*this);
	blas::copt_blas_axpy(this->size(),-1.0,mat.dataPtr(),1,result.dataPtr(),1);
	return result;
}

template<class scalar,class index>
VectorBase<scalar,index> MatrixBase<scalar,index>::operator*(const VectorBase<scalar,index>& vec )const
{
	if ( __cols != vec.size() )
		throw COException("MatrixBase multiply error: the size of MatrixBase and vector are not consistent!");
	VectorBase<scalar,index> result(__rows);
	blas::copt_blas_gemv(CblasColMajor,CblasNoTrans,__rows,__cols,1.0,this->dataPtr(),__rows,vec.dataPtr(),1,0.0,result.dataPtr(),1);
	return result;
}

template<class scalar,class index>
MatrixBase<scalar,index> MatrixBase<scalar,index>::operator*(const MatrixBase& mat )const
{
	if ( __cols != mat.rows() )
		throw COException("MatrixBase multiply error: the size of two matrices are not consistent!");
	MatrixBase result(__rows,mat.cols());
	blas::copt_blas_gemm(CblasColMajor,CblasNoTrans,CblasNoTrans,__rows,mat.cols(),__cols,1.0,this->dataPtr(),__rows,mat.dataPtr(),__cols,0.0,result.dataPtr(),__rows);
	return result;
}

template<class scalar,class index>
MatrixBase<scalar,index> MatrixBase<scalar,index>::identity(
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


template<class scalar,class index>
void MatrixBase<scalar,index>::blockFromMatrix(const MatrixBase& mat,const std::set<index>& rownums,const std::set<index>& colnums)
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

template<class scalar,class index>
void MatrixBase<scalar,index>::columnBlockFromMatrix(const MatrixBase& mat,const std::set<index>& colnums)
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

template<class scalar,class index>
void MatrixBase<scalar,index>::rowBlockFromMatrix(const MatrixBase& mat,const std::set<index>& rownums)
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

template<class scalar,class index>
void MatrixBase<scalar,index>::blockFromMatrix(const MatrixBase& mat,const std::vector<index>& rownums,const std::vector<index>& colnums)
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

template<class scalar,class index> template<class InputIterator>
void MatrixBase<scalar,index>::blockFromMatrix(
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

template<class scalar,class index>
void MatrixBase<scalar,index>::columnBlockFromMatrix(const MatrixBase& mat,const std::vector<index>& colnums)
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

template<class scalar,class index> template<class InputIterator>
void MatrixBase<scalar,index>::columnBlockFromMatrix(
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

template<class scalar,class index>
void MatrixBase<scalar,index>::rowBlockFromMatrix(const MatrixBase& mat,const std::vector<index>& rownums)
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
void MatrixBase<scalar,index>::rowBlockFromMatrix(
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

template<class scalar,class index>
void MatrixBase<scalar,index>::combineAlongRow(const MatrixBase& m1,const MatrixBase& m2)
{
	MatrixBase::stCombineAlongRow(m1,m2,*this);
}

template<class scalar,class index>
void MatrixBase<scalar,index>::combineAlongColumn(const MatrixBase& m1,const MatrixBase& m2)
{
	MatrixBase::stCombineAlongColumn(m1,m2,*this);
}

template<class scalar,class index>
void MatrixBase<scalar,index>::stCombineAlongRow(const MatrixBase& m1,const MatrixBase& m2,MatrixBase& m)
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

template<class scalar,class index>
void MatrixBase<scalar,index>::stCombineAlongColumn(const MatrixBase& m1,const MatrixBase& m2,MatrixBase& m)
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

template<class scalar,class index>
void MatrixBase<scalar,index>::setRandom(const index rows,const index cols)
{
	if(rows<0||cols<0)
		throw COException("Please make sure that the number of row and column is bigger than zero!");
	std::mt19937 eng(time(NULL));
	std::uniform_real_distribution<scalar> unif(0.0,1.0);
	this->resize(rows,cols);
	for ( int i = 0 ; i < rows ; ++ i )
		for ( int j = 0 ; j < cols ; ++ j )
			this->operator()(i,j)=unif(eng);
}

template<class scalar,class index>
MatrixBase<scalar,index> MatrixBase<scalar,index>::random(const index rows,const index cols)
{
	MatrixBase result;
	result.setRandom(rows,cols);
	return result;
}

template<class scalar,class index>
void MatrixBase<scalar,index>::mtm( MatrixBase& mat ) const
{
	int m = this->rows();
	int n = this->cols();
	mat.resize(n,n);
	blas::copt_blas_syrk(CblasColMajor,CblasUpper,CblasTrans,n,m,1.0,this->dataPtr(),m,0.0,mat.dataPtr(),n);
	mat.setSymmetricFlag(true);
}

/*******************Implementation of Triplet******************/

template<class scalar,class index>
TripletBase<scalar,index>::TripletBase(
	const index r,
	const index c,
	const scalar v)
	:
	__r(r),
	__c(c),
	__v(v)
{
}

template<class scalar,class index>
TripletBase<scalar,index>::~TripletBase()
{
}

template<class scalar,class index>
const index& TripletBase<scalar,index>::rowIndex() const
{
	return __r;
}

template<class scalar,class index>
const index& TripletBase<scalar,index>::columnIndex() const
{
	return __c;
}

template<class scalar,class index>
const scalar& TripletBase<scalar,index>::value() const
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

template<class scalar,class index>
void SpMatrixBase<scalar,index>::judgeRationality()
{
	for ( index i = 0 ; i < __elesize ; ++ i )
	{
		if ( __rowind[i] >= __rows )
			throw COException("Sparse matrix not rational: row index out of range!");
	}
}

template<class scalar,class index>
const typename SpMatrixBase<scalar,index>::scalar SpMatrixBase<scalar,index>::__zero = static_cast<scalar>(0.0);

template<class scalar,class index>
SpMatrixBase<scalar,index>::SpMatrixBase()
	:
	__rows(0),
	__cols(0),
	__elesize(0),
	__colptr(NULL),
	__rowind(NULL),
	__vals(NULL)
{
}

template<class scalar,class index>
SpMatrixBase<scalar,index>::SpMatrixBase(
	const index 						rows,
	const index 						cols,
	const index 						elesize,
	const index*						colptr,
	const index*						rowind,
	const scalar*						vals)
	:
	__rows(rows),
	__cols(cols),
	__elesize(elesize),
	__colptr(new index[cols+1]),
	__rowind(new index[elesize]),
	__vals(new scalar[elesize])
{

	blas::copt_blas_copy(__cols+1,colptr,1,__colptr,1);
	blas::copt_blas_copy(__elesize,rowind,1,__rowind,1);
	blas::copt_blas_copy(__elesize,vals,1,__vals,1);

	judgeRationality();
}

template<class scalar,class index>
SpMatrixBase<scalar,index>::SpMatrixBase(
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

template<class scalar,class index>
SpMatrixBase<scalar,index>::~SpMatrixBase()
{
	if(__rowind)
		SAFE_DELETE_ARRAY(__rowind);
	if(__vals)
		SAFE_DELETE_ARRAY(__vals);
	if(__colptr)
		SAFE_DELETE_ARRAY(__colptr);
}

template<class scalar,class index>
void SpMatrixBase<scalar,index>::setSparseMatrix(
	const index 					rows,
	const index 					cols,
	const index 					elesize,
	const index*			 		colptr,
	const index*					rowind,
	const scalar*			 		vals)
{
	clear();
	__rows = rows;
	__cols = cols;
	__elesize = elesize;
	
	__rowind = new index[__elesize];
	__vals = new scalar[__elesize];
	__colptr = new index[__cols+1];

	blas::copt_blas_copy(__elesize,rowind,1,__rowind,1);
	blas::copt_blas_copy(__elesize,vals,1,__vals,1);
	blas::copt_blas_copy(__cols+1,colptr,1,__colptr,1);

	judgeRationality();
}

template<class scalar,class index>
SpMatrixBase<scalar,index>& SpMatrixBase<scalar,index>::operator=(const SpMatrixBase& mat)
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

template<class scalar,class index>
void SpMatrixBase<scalar,index>::setFromTriplets(
	const index rows,
	const index cols,
	std::vector<Triplet>& triplets)
{
	clear();
	__rows = rows;
	__cols = cols;
	// sort the triplets according to the column index at first
	std::sort(triplets.begin(),triplets.end(),columnComparison<Triplet>());
	// compute how many elements there are in one column
	index colind = 0 , ip = 0;
	for ( index i = 0 ; i < triplets.size() ; ++ i )
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
	std::vector<index> 	colcounts(__cols,0);
	std::list<index> 		rowinds;
	std::list<scalar>	vals;
	colind = 0 , ip = 0;
	index count = 0;
	for ( index i = 0 ; i < triplets.size() ; ++ i )
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

	__colptr = new index[__cols+1];
	index columncount = 0;
	for ( index i = 0 ; i < __cols ; ++ i )
	{
		__colptr[i] = columncount;
		columncount += colcounts[i];
	}
	__colptr[__cols] = columncount;

	__rowind = new index[rowinds.size()];
	index i = 0;
	for ( typename std::list<index>::iterator iter = rowinds.begin() ; iter != rowinds.end() ; ++ iter , ++ i )
		__rowind[i] = *iter;

	__vals = new scalar[vals.size()];
	i = 0;
	for ( typename std::list<scalar>::iterator iter = vals.begin() ; iter != vals.end() ; ++ iter , ++ i )
		__vals[i] = *iter;

	__elesize = rowinds.size();
}

template<class scalar,class index> template<class InputIterator>
void SpMatrixBase<scalar,index>::setFromTriplets(
	const index rows,
	const index cols,
	const InputIterator& begin,
	const InputIterator& end)
{
	clear();
	__rows = rows;
	__cols = cols;
	// sort the triplets according to the column index at first
	std::sort(begin,end,columnComparison<Triplet>());
	// compute how many elements there are in one column
	index colind = 0 , ip = 0;
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
	std::vector<index> 		colcounts(__cols,0);
	std::list<index> 		rowinds;
	std::list<scalar>	vals;
	colind = 0 , ip = 0;
	previter = begin;
	index count = 0;
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
		else if(iter->rowIndex()==previter->rowIndex()&&iter->columnIndex()==previter->columnIndex())
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

	__colptr = new index[__cols+1];
	index columncount = 0;
	for ( index i = 0 ; i < __cols ; ++ i )
	{
		__colptr[i] = columncount;
		columncount += colcounts[i];
	}
	__colptr[__cols] = columncount;

	__rowind = new index[rowinds.size()];
	index i = 0;
	for ( typename std::list<index>::iterator iter = rowinds.begin() ; iter != rowinds.end() ; ++ iter , ++ i )
		__rowind[i] = *iter;

	__vals = new scalar[vals.size()];
	i = 0;
	for ( typename std::list<scalar>::iterator iter = vals.begin() ; iter != vals.end() ; ++ iter , ++ i )
		__vals[i] = *iter;

	__elesize = rowinds.size();
}

template<class scalar,class index> template<class InputIterator>
void SpMatrixBase<scalar,index>::fastSetFromTriplets(
	const index rows,
	const index cols,
	const InputIterator& begin,
	const InputIterator& end)
{
	clear();
	__rows = rows;
	__cols = cols;
	__colptr = new index[__cols+1];
	// compute nnz
	index nnz = 0;
	for ( InputIterator iter = begin ; iter != end ; ++ iter )
		++ nnz;
	__elesize = nnz;
	__rowind = new index[nnz];
	__vals = new scalar[nnz];
	InputIterator previter = begin;
	index count = 0 , i = 0;
	std::vector<index> counts (__cols,0);
	for ( InputIterator iter = begin ; iter != end ; ++ iter , ++ i )
	{
		if(iter->columnIndex()!=previter->columnIndex())
		{
			counts[previter->columnIndex()] = count;
			count = 0;
			previter = iter;
		}
		__rowind[i] = iter->rowIndex();
		__vals[i] = iter->value();
		++ count;
	}
	counts[previter->columnIndex()] = count;
	count = 0;
	for ( index c = 0 ; c < __cols ; ++ c )
	{
		__colptr[c] = count;
		count += counts[c];
	}
	__colptr[__cols] = count;
}

template<class scalar,class index>
void SpMatrixBase<scalar,index>::clear()
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

template<class scalar,class index>
const index& SpMatrixBase<scalar,index>::rows() const
{
	return __rows;
}

template<class scalar,class index>
const index& SpMatrixBase<scalar,index>::cols() const
{
	return __cols;
}

template<class scalar,class index>
const index& SpMatrixBase<scalar,index>::elementSize() const
{
	return __elesize;
}

template<class scalar,class index>
const index* SpMatrixBase<scalar,index>::columnPointer() const
{
	return __colptr;
}

template<class scalar,class index>
const index* SpMatrixBase<scalar,index>::rowIndex() const
{
	return __rowind;
}

template<class scalar,class index>
const scalar* SpMatrixBase<scalar,index>::values() const
{
	return __vals;
}

template<class scalar,class index>
void SpMatrixBase<scalar,index>::scale(const scalar s)
{
	for ( index i = 0 ; i < __elesize ; ++ i )
		__vals[i] *= s;
}

template<class scalar,class index>
void SpMatrixBase<scalar,index>::neg()
{
	for ( index i = 0 ; i < __elesize ; ++ i )
		__vals[i] = -__vals[i];
}

template<class scalar,class index>
const scalar& SpMatrixBase<scalar,index>::operator()(
	const index i,
	const index j) const
{
	if(i>=__rows||j>=__cols)
		throw COException("Sparse Matrix error, index out of range!");
	index ip = __colptr[j],in=__colptr[j+1];
	for ( index ind = ip ; ind < in ; ++ ind )
	{
		if( i == __rowind[ind] )
			return __vals[ind];
	}
	return __zero;
}

template<class scalar,class index>
SpMatrixBase<scalar,index> SpMatrixBase<scalar,index>::operator*(const SpMatrixBase& mat) const
{
	if(__cols != mat.rows() )
		throw COException("Multiplication error: matrix index does not fit!");
	std::list<Triplet> tris;
	for ( index c = 0 ; c < __cols ; ++ c )
	{
		index ci = __colptr[c] , cn = __colptr[c+1];
		for ( index r = ci ; r < cn ; ++ r )
		{
			index rind = __rowind[r];
			// traverse mat
			for ( index mc = 0 ; mc < mat.cols() ; ++ mc )
			{
				index mci = mat.columnPointer()[mc], mcn = mat.columnPointer()[mc+1];
				for ( index mr = mci ; mr < mcn ; ++ mr )
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


template<class scalar,class index>
SpMatrixBase<scalar,index> SpMatrixBase<scalar,index>::operator+ ( const SpMatrixBase& mat ) const
{
	if( __cols != mat.cols() || __rows != mat.rows() )
		throw COException("Sparse matrix summation error: index does not fit!");
	std::list<index> inds;
	std::list<scalar> vals;
	index *colptr = new index[__cols+1];
	index count = 0;
	for ( index c = 0 ; c < __cols ; ++ c )
	{
		colptr[c] = count;
		index ci1 = __colptr[c],cn1 = __colptr[c+1],ci2 = mat.columnPointer()[c],cn2=mat.columnPointer()[c+1];
		index i1 = ci1,i2=ci2;
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
	index *rowind = new index[inds.size()];
	index i = 0;
	for ( typename std::list<index>::iterator iter = inds.begin() ; iter != inds.end() ; ++ iter , ++i )
		rowind[i] = *iter;
	i = 0;
	scalar *vs = new scalar[vals.size()];
	for ( typename std::list<scalar>::iterator iter = vals.begin() ; iter != vals.end() ; ++ iter , ++ i )
		vs[i] = *iter;
	SpMatrixBase result;
	result.setSparseMatrix(__rows,__cols,inds.size(),colptr,rowind,vs);

	delete[]rowind;
	delete[]colptr;
	delete[]vs;
	return result;
}

template<class scalar,class index>
SpMatrixBase<scalar,index> SpMatrixBase<scalar,index>::operator- ( const SpMatrixBase& mat ) const
{
	if( __cols != mat.cols() || __rows != mat.rows() )
		throw COException("Sparse matrix summation error: index does not fit!");
	std::list<index> inds;
	std::list<scalar> vals;
	index *colptr = new index[__cols+1];
	index count = 0;
	for ( index c = 0 ; c < __cols ; ++ c )
	{
		colptr[c] = count;
		index ci1 = __colptr[c],cn1 = __colptr[c+1],ci2 = mat.columnPointer()[c],cn2=mat.columnPointer()[c+1];
		index i1 = ci1,i2=ci2;
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
	index *rowind = new index[inds.size()];
	index i = 0;
	for ( typename std::list<index>::iterator iter = inds.begin() ; iter != inds.end() ; ++ iter , ++i )
		rowind[i] = *iter;
	i = 0;
	scalar *vs = new scalar[vals.size()];
	for ( typename std::list<scalar>::iterator iter = vals.begin() ; iter != vals.end() ; ++ iter , ++ i )
		vs[i] = *iter;
	SpMatrixBase result;
	result.setSparseMatrix(__rows,__cols,inds.size(),colptr,rowind,vs);
	delete[]rowind;
	delete[]colptr;
	delete[]vs;
	return result;
}

template<class scalar,class index>
SpMatrixBase<scalar,index> SpMatrixBase<scalar,index>::operator- () const
{
	SpMatrixBase result(*this);
	result.neg();
	return result;
}

template<class scalar,class index>
typename SpMatrixBase<scalar,index>::Vector SpMatrixBase<scalar,index>::operator*(const Vector& vec) const
{
	if(__cols!=vec.size() )
		throw COException("Multiplication error: index does not fit!");
	Vector result(__rows);
	for ( index i = 0 ; i < __cols ; ++ i )
	{
		index ip = __colptr[i] , in = __colptr[i+1];
		for ( index r = ip ; r < in ; ++ r ){
			result[__rowind[r]] += __vals[r]*vec[i];
		}
	}
	return result;
}


template<class scalar,class index, class T>
SpMatrixBase<scalar,index> operator* (const T s,const SpMatrixBase<scalar,index>& mat)
{
	return mat.operator*(static_cast<scalar>(s));
}

template<class scalar,class index>
MatrixBase<scalar,index> SpMatrixBase<scalar,index>::toDenseMatrix() const
{
	MatrixBase<scalar,index> result(__rows,__cols);
	for ( index i = 0 ; i < __cols ; ++ i )
	{
		index ip = __colptr[i] , in = __colptr[i+1];
		for ( index r = ip ; r < in ; ++ r )
		{
			result(__rowind[r],i) = __vals[r];
		}
	}
	return result;
}

template<class scalar,class index>
VectorBase<scalar,index> SpMatrixBase<scalar,index>::solve(const VectorBase<scalar,index>& vec)
{
	return UMFLinearSolver<SpMatrixBase>(*this).solve(vec);
}

template<class scalar,class index>
SpMatrixBase<scalar,index> SpMatrixBase<scalar,index>::operator* ( const scalar s ) const
{
	SpMatrixBase result(*this);
	result.scale(s);
	return result;
}

template<class scalar,class index>
SpMatrixBase<scalar,index> operator*(const scalar s,const SpMatrixBase<scalar,index>& mat)
{
	return mat.operator*(s);
}

}// End of namespace COPT

#endif