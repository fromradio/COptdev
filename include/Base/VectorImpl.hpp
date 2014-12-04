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
template<class scalar,class index>
VectorBase<scalar,index>::VectorBase()
	:Array<scalar,index>()
{
}

template<class scalar,class index>
scalar& VectorBase<scalar,index>::operator() ( const index i )
{
	return this->operator[](i);
}

template<class scalar,class index>
const scalar& VectorBase<scalar,index>::operator() (const index i ) const
{
	return this->operator[](i);
}

template<class scalar,class index>
bool VectorBase<scalar,index>::operator<(const VectorBase<scalar,index>& vec)const
{
	if (this->size() != vec.size() )
		return false;
	for ( index i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)>=vec[i])
			return false;
	}
	return true;
}

template<class scalar,class index>
bool VectorBase<scalar,index>::operator<=(const VectorBase<scalar,index>& vec)const
{
	if (this->size() != vec.size() )
		return false;
	for ( index i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)>vec[i])
			return false;
	}
	return true;
}

template<class scalar,class index>
bool VectorBase<scalar,index>::operator>(const VectorBase<scalar,index>& vec)const
{
	if (this->size() != vec.size() )
		return false;
	for ( index i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)<=vec[i])
			return false;
	}
	return true;
}

template<class scalar,class index>
bool VectorBase<scalar,index>::operator>=(const VectorBase<scalar,index>& vec)const
{
	if (this->size() != vec.size() )
		return false;
	for ( index i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)<vec[i])
			return false;
	}
	return true;
}

template<class scalar,class index>
bool VectorBase<scalar,index>::operator==(const VectorBase<scalar,index>& vec)const
{
	if (this->size() != vec.size() )
		return false;
	for ( index i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)!=vec[i])
			return false;
	}
	return true;
}

template<class scalar,class index>
bool VectorBase<scalar,index>::operator!=(const VectorBase<scalar,index>& vec)const
{
	if ( this->size() != vec.size() )
		return true;
	for ( index i = 0 ; i < this->size() ; ++ i ){
		if(this->operator[](i)!=vec[i])
			return true;
	}
	return false;
}

template<class scalar,class index>
scalar VectorBase<scalar,index>::normalize()
{
	scalar norm = std::sqrt(squaredNorm());
	scale(1.0/norm);
	return norm;
}

template<class scalar,class index>
MatrixBase<scalar,index> VectorBase<scalar,index>::mulTrans(const VectorBase<scalar,index>& vec) const
{
	index 					m = this->size();
	index 					n = vec.size();
	MatrixBase<scalar,index> 	result(m,n);
	for ( index i = 0 ; i < m ; ++ i )
		for ( index j = 0 ; j < n ; ++ j )
			result(i,j) = this->operator[](i)*vec[j];
	return result;
}


/*			the transpose of a matrix multiplies a vector
 */
template<class scalar,class index>
VectorBase<scalar,index> VectorBase<scalar,index>::transMul(const MatrixBase<scalar,index>& mat) const
{
	return mat.transpose()*(*this);
}

template<class scalar,class index>
VectorBase<scalar,index> VectorBase<scalar,index>::block(const std::set<index>& indices) const
{
	if( *indices.rbegin() >= this->size() )
	{
		throw COException("Index out of range in Vector blocking!");
	}
	VectorBase<scalar,index> result(indices.size());
	index i = 0;
	for ( typename std::set<index>::const_iterator iter = indices.begin() ; iter != indices.end() ; ++ iter ){
		result[i] = this->operator[](*iter);
		++ i;
	} 
	return result;
}

template<class scalar,class index>
void VectorBase<scalar,index>::blockFromVector(const VectorBase& vec,const std::set<index>& indices)
{
	if( *indices.rbegin() >= vec.size() )
	{
		throw COException("Index out of range in Vector blocking!");
	}
	this->resize(indices.size());
	index i = 0;
	for ( typename std::set<index>::const_iterator iter = indices.begin() ; iter != indices.end() ; ++ iter ){
		this->operator[](i) = vec[*iter];
		++ i;
	}
}

template<class scalar,class index>
VectorBase<scalar,index> VectorBase<scalar,index>::block(const std::vector<index>& indices) const
{
	VectorBase<scalar,index> result(indices.size());
	for ( index i = 0 ; i < indices.size() ; ++ i ){
		if (indices[i] >= this->size())
		{
			throw COException("Index out of range in Vector blocking!");
		}
		result[i] = this->operator[](indices[i]);
	}
	return result;
}

template<class scalar,class index>
void VectorBase<scalar,index>::blockFromVector(const VectorBase& vec,const std::vector<index>& indices)
{
	this->resize(indices.size());
	for ( index i = 0 ; i < indices.size() ; ++ i ){
		if(indices[i] >= vec.size() )
		{
			throw COException("Index out of range in Vector blocking!");
		}
		this->operator[](i)=vec[indices[i]];
	}
}

template<class scalar,class index>
void VectorBase<scalar,index>::combine(const VectorBase& v1,const VectorBase& v2)
{
	VectorBase::stCombine(v1,v2,*this);
}

template<class scalar,class index>
void VectorBase<scalar,index>::stCombine(const VectorBase& v1,const VectorBase& v2,VectorBase& v)
{
	v.resize(v1.size()+v2.size());
	index n = v1.size();
	for ( index i = 0 ; i < n ; ++ i )
		v[i] = v1[i];
	for ( index i = 0 ; i < v2.size() ; ++ i )
		v[i+n] = v2[i];
}

}// End of namespace COPT


#endif