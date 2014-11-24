// 		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef VECTOR_IMPL_H
#define VECTOR_IMPL_H
/*
				Implementation of class 'VectorBase' taking number type as template
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

template<class ScalarType>
ScalarType& VectorBase<ScalarType>::operator() ( const size_t i )
{
	return this->operator[](i);
}

template<class ScalarType>
const ScalarType& VectorBase<ScalarType>::operator() (const size_t i ) const
{
	return this->operator[](i);
}

template<class ScalarType>
bool VectorBase<ScalarType>::operator<(const VectorBase<ScalarType>& vec)const
{
	if (this->size() != vec.size() )
		return false;
	for ( int i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)>=vec[i])
			return false;
	}
	return true;
}

template<class ScalarType>
bool VectorBase<ScalarType>::operator<=(const VectorBase<ScalarType>& vec)const
{
	if (this->size() != vec.size() )
		return false;
	for ( int i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)>vec[i])
			return false;
	}
	return true;
}

template<class ScalarType>
bool VectorBase<ScalarType>::operator>(const VectorBase<ScalarType>& vec)const
{
	if (this->size() != vec.size() )
		return false;
	for ( int i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)<=vec[i])
			return false;
	}
	return true;
}

template<class ScalarType>
bool VectorBase<ScalarType>::operator>=(const VectorBase<ScalarType>& vec)const
{
	if (this->size() != vec.size() )
		return false;
	for ( int i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)<vec[i])
			return false;
	}
	return true;
}

template<class ScalarType>
bool VectorBase<ScalarType>::operator==(const VectorBase<ScalarType>& vec)const
{
	if (this->size() != vec.size() )
		return false;
	for ( int i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)!=vec[i])
			return false;
	}
	return true;
}

template<class ScalarType>
bool VectorBase<ScalarType>::operator!=(const VectorBase<ScalarType>& vec)const
{
	if ( this->size() != vec.size() )
		return true;
	for ( int i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)!=vec[i])
			return true;
	}
	return false;
}

template<class ScalarType>
MatrixBase<ScalarType> VectorBase<ScalarType>::mulTrans(const VectorBase<ScalarType>& vec) const
{
	size_t 					m = this->size();
	size_t 					n = vec.size();
	MatrixBase<ScalarType> 	result(m,n);
	for ( int i = 0 ; i < m ; ++ i )
		for ( int j = 0 ; j < n ; ++ j )
			result(i,j) = this->operator[](i)*vec[j];
	return result;
}


/*			the transpose of a matrix multiplies a vector
 */
template<class ScalarType>
VectorBase<ScalarType> VectorBase<ScalarType>::transMul(const MatrixBase<ScalarType>& mat) const
{
	return mat.transpose()*(*this);
}

template<class ScalarType>
VectorBase<ScalarType> VectorBase<ScalarType>::block(const std::set<size_t>& indices) const
{
	if( *indices.rbegin() >= this->size() )
	{
		throw COException("Index out of range in Vector blocking!");
	}
	VectorBase<ScalarType> result(indices.size());
	int i = 0;
	for ( std::set<size_t>::const_iterator iter = indices.begin() ; iter != indices.end() ; ++ iter ){
		result[i] = this->operator[](*iter);
		++ i;
	} 
	return result;
}

template<class ScalarType>
void VectorBase<ScalarType>::blockFromVector(const VectorBase& vec,const std::set<size_t>& indices)
{
	if( *indices.rbegin() >= vec.size() )
	{
		throw COException("Index out of range in Vector blocking!");
	}
	this->resize(indices.size());
	int i = 0;
	for ( std::set<size_t>::const_iterator iter = indices.begin() ; iter != indices.end() ; ++ iter ){
		this->operator[](i) = vec[*iter];
		++ i;
	}
}

template<class ScalarType>
VectorBase<ScalarType> VectorBase<ScalarType>::block(const std::vector<size_t>& indices) const
{
	VectorBase<ScalarType> result(indices.size());
	for ( int i = 0 ; i < indices.size() ; ++ i ){
		if (indices[i] >= this->size())
		{
			throw COException("Index out of range in Vector blocking!");
		}
		result[i] = this->operator[](indices[i]);
	}
	return result;
}

template<class ScalarType>
void VectorBase<ScalarType>::blockFromVector(const VectorBase& vec,const std::vector<size_t>& indices)
{
	this->resize(indices.size());
	for ( int i = 0 ; i < indices.size() ; ++ i ){
		if(indices[i] >= vec.size() )
		{
			throw COException("Index out of range in Vector blocking!");
		}
		this->operator[](i)=vec[indices[i]];
	}
}

template<class ScalarType>
void VectorBase<ScalarType>::combine(const VectorBase& v1,const VectorBase& v2)
{
	VectorBase::stCombine(v1,v2,*this);
}

template<class ScalarType>
void VectorBase<ScalarType>::stCombine(const VectorBase& v1,const VectorBase& v2,VectorBase& v)
{
	v.resize(v1.size()+v2.size());
	size_t n = v1.size();
	for ( size_t i = 0 ; i < n ; ++ i )
		v[i] = v1[i];
	for ( size_t i = 0 ; i < v2.size() ; ++ i )
		v[i+n] = v2[i];
}

}// End of namespace COPT


#endif