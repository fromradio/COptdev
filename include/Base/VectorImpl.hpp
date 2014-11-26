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
template<class ScalarType,class Size>
VectorBase<ScalarType,Size>::VectorBase()
	:Array<ScalarType>()
{
}

template<class ScalarType,class Size>
ScalarType& VectorBase<ScalarType,Size>::operator() ( const Size i )
{
	return this->operator[](i);
}

template<class ScalarType,class Size>
const ScalarType& VectorBase<ScalarType,Size>::operator() (const Size i ) const
{
	return this->operator[](i);
}

template<class ScalarType,class Size>
bool VectorBase<ScalarType,Size>::operator<(const VectorBase<ScalarType,Size>& vec)const
{
	if (this->size() != vec.size() )
		return false;
	for ( int i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)>=vec[i])
			return false;
	}
	return true;
}

template<class ScalarType,class Size>
bool VectorBase<ScalarType,Size>::operator<=(const VectorBase<ScalarType,Size>& vec)const
{
	if (this->size() != vec.size() )
		return false;
	for ( int i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)>vec[i])
			return false;
	}
	return true;
}

template<class ScalarType,class Size>
bool VectorBase<ScalarType,Size>::operator>(const VectorBase<ScalarType,Size>& vec)const
{
	if (this->size() != vec.size() )
		return false;
	for ( int i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)<=vec[i])
			return false;
	}
	return true;
}

template<class ScalarType,class Size>
bool VectorBase<ScalarType,Size>::operator>=(const VectorBase<ScalarType,Size>& vec)const
{
	if (this->size() != vec.size() )
		return false;
	for ( int i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)<vec[i])
			return false;
	}
	return true;
}

template<class ScalarType,class Size>
bool VectorBase<ScalarType,Size>::operator==(const VectorBase<ScalarType,Size>& vec)const
{
	if (this->size() != vec.size() )
		return false;
	for ( int i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)!=vec[i])
			return false;
	}
	return true;
}

template<class ScalarType,class Size>
bool VectorBase<ScalarType,Size>::operator!=(const VectorBase<ScalarType,Size>& vec)const
{
	if ( this->size() != vec.size() )
		return true;
	for ( int i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)!=vec[i])
			return true;
	}
	return false;
}

template<class ScalarType,class Size>
MatrixBase<ScalarType,Size> VectorBase<ScalarType,Size>::mulTrans(const VectorBase<ScalarType,Size>& vec) const
{
	Size 					m = this->size();
	Size 					n = vec.size();
	MatrixBase<ScalarType,Size> 	result(m,n);
	for ( int i = 0 ; i < m ; ++ i )
		for ( int j = 0 ; j < n ; ++ j )
			result(i,j) = this->operator[](i)*vec[j];
	return result;
}


/*			the transpose of a matrix multiplies a vector
 */
template<class ScalarType,class Size>
VectorBase<ScalarType,Size> VectorBase<ScalarType,Size>::transMul(const MatrixBase<ScalarType,Size>& mat) const
{
	return mat.transpose()*(*this);
}

template<class ScalarType,class Size>
VectorBase<ScalarType,Size> VectorBase<ScalarType,Size>::block(const std::set<Size>& indices) const
{
	if( *indices.rbegin() >= this->size() )
	{
		throw COException("Index out of range in Vector blocking!");
	}
	VectorBase<ScalarType,Size> result(indices.size());
	int i = 0;
	for ( typename std::set<Size>::const_iterator iter = indices.begin() ; iter != indices.end() ; ++ iter ){
		result[i] = this->operator[](*iter);
		++ i;
	} 
	return result;
}

template<class ScalarType,class Size>
void VectorBase<ScalarType,Size>::blockFromVector(const VectorBase& vec,const std::set<Size>& indices)
{
	if( *indices.rbegin() >= vec.size() )
	{
		throw COException("Index out of range in Vector blocking!");
	}
	this->resize(indices.size());
	int i = 0;
	for ( typename std::set<Size>::const_iterator iter = indices.begin() ; iter != indices.end() ; ++ iter ){
		this->operator[](i) = vec[*iter];
		++ i;
	}
}

template<class ScalarType,class Size>
VectorBase<ScalarType,Size> VectorBase<ScalarType,Size>::block(const std::vector<Size>& indices) const
{
	VectorBase<ScalarType,Size> result(indices.size());
	for ( int i = 0 ; i < indices.size() ; ++ i ){
		if (indices[i] >= this->size())
		{
			throw COException("Index out of range in Vector blocking!");
		}
		result[i] = this->operator[](indices[i]);
	}
	return result;
}

template<class ScalarType,class Size>
void VectorBase<ScalarType,Size>::blockFromVector(const VectorBase& vec,const std::vector<Size>& indices)
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

template<class ScalarType,class Size>
void VectorBase<ScalarType,Size>::combine(const VectorBase& v1,const VectorBase& v2)
{
	VectorBase::stCombine(v1,v2,*this);
}

template<class ScalarType,class Size>
void VectorBase<ScalarType,Size>::stCombine(const VectorBase& v1,const VectorBase& v2,VectorBase& v)
{
	v.resize(v1.size()+v2.size());
	Size n = v1.size();
	for ( Size i = 0 ; i < n ; ++ i )
		v[i] = v1[i];
	for ( Size i = 0 ; i < v2.size() ; ++ i )
		v[i+n] = v2[i];
}

}// End of namespace COPT


#endif