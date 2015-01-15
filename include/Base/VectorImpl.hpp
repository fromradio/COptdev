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


#ifndef VECTOR_IMPL_HPP__
#define VECTOR_IMPL_HPP__

/********************Implementation of class 'VectorBase'******************/
namespace COPT
{

template<class scalar,class index>
VectorBase<scalar,index>::VectorBase()
	:
	Array<scalar,index>()
{
}

template<class scalar,class index>
VectorBase<scalar,index>::VectorBase( const index size , scalar *data )
	:
	Arr(size,data)
{
}

template<class scalar,class index>
VectorBase<scalar,index>::VectorBase( const VectorBase& vec )
	:
	Arr()
{
	if (vec.isReferred())
	{
		// the vector is a referred vector
		this->setReferredArray(vec.size(),const_cast<scalar*>(vec.dataPtr()),vec.interval());
	}
	else{
		// copy the vector data
		this->setArray(vec.size(),vec.dataPtr(),vec.interval());
	}
}

template<class scalar,class index>
VectorBase<scalar,index>::VectorBase( const index size, const referred_array& tag, scalar *data, const index inter )
	:
	Arr(size,tag,data,inter)
{
}

template<class scalar,class index>
VectorBase<scalar,index>::VectorBase( const index size , const referred_array& tag , const scalar* data ,const index inter )
	:
	Arr(size,tag,const_cast<scalar*>(data),inter)
{
}

template<class scalar,class index>
VectorBase<scalar,index>::VectorBase(const std::vector<scalar>& vec)
	:
	Arr()
{
	this->resize(vec.size(),1);
	for ( index i = 0 ; i < this->size() ; ++ i )
		this->operator[](i) = vec[i];
}

template<class scalar,class index>
VectorBase<scalar,index>::~VectorBase()
{
}

template<class scalar,class index>
VectorBase<scalar,index>& VectorBase<scalar,index>::operator= ( const VectorBase& vec )
{
	if (vec.isReferred()){
		this->setReferredArray(vec.size(),const_cast<scalar*>(vec.dataPtr()),vec.interval());
	}
	else{
		this->setArray(vec.size(),vec.dataPtr(),vec.interval());
	}
	return *this;
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
		if(LargerThan(this->operator[](i),vec[i]))
			return false;
	}
	return true;
}

template<class scalar,class index>
bool VectorBase<scalar,index>::operator<(const scalar s)const
{
	for ( index i = 0 ; i < this->size() ; ++ i )
	{
		if(LargerThan(this->operator[](i),s))
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
		if(StrictLargerThan(this->operator[](i),vec[i]))
			return false;
	}
	return true;
}

template<class scalar,class index>
bool VectorBase<scalar,index>::operator<=(const scalar s)const
{
	for ( index i = 0 ; i < this->size() ; ++ i )
	{
		if(StrictLargerThan(this->operator[](i),s))
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
		if(LessThan(this->operator[](i),vec[i]))
			return false;
	}
	return true;
}

template<class scalar,class index>
bool VectorBase<scalar,index>::operator>(const scalar s)const
{
	for ( index i = 0 ; i < this->size() ; ++ i )
	{
		if(LessThan(this->operator[](i),s))
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
		if(StrictLessThan(this->operator[](i),vec[i]))
			return false;
	}
	return true;
}

template<class scalar,class index>
bool VectorBase<scalar,index>::operator>=(const scalar s)const
{
	for ( index i = 0 ; i < this->size() ; ++ i )
	{
		if(StrictLessThan(this->operator[](i),s))
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
		if(!IS_ZERO(std::abs(this->operator[](i)-vec[i])))
			return false;
	}
	return true;
}

template<class scalar,class index>
bool VectorBase<scalar,index>::operator==(const scalar s)const
{
	for ( index i = 0 ; i < this->size() ; ++ i )
	{
		if(!IS_ZERO(std::abs(this->operator[](i)-s)))
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
bool VectorBase<scalar,index>::operator!=(const scalar s)const
{
	for ( index i = 0 ; i < this->size() ; ++ i )
	{
		if(this->operator[](i)!=s)
			return true;
	}
	return false;
}

template<class scalar,class index>
typename VectorBase<scalar,index>::podscalar VectorBase<scalar,index>::squaredNorm() const
{
	podscalar norm = blas::copt_blas_nrm2(this->size(),this->dataPtr(),this->interval());
	return norm*norm;
}

template<class scalar,class index>
auto VectorBase<scalar,index>::norm()const->podscalar
{
	podscalar norm = blas::copt_blas_nrm2(this->size(),this->dataPtr(),this->interval());
	return norm;
}

template<class scalar,class index>
typename VectorBase<scalar,index>::podscalar VectorBase<scalar,index>::absNorm() const
{
	podscalar result = 0;
	for (index i = 0 ; i < this->size() ; ++ i )
	{
		result += std::abs(this->operator[](i));
	}
	return result;
}

template<class scalar,class index>
typename VectorBase<scalar,index>::podscalar VectorBase<scalar,index>::normalize()
{
	podscalar norm = std::sqrt(squaredNorm());
	scale(1.0/norm);
	return norm;
}

template<class scalar,class index>
typename VectorBase<scalar,index>::scalar VectorBase<scalar,index>::dot(const VectorBase& vec) const
{
	if(this->size()!=vec.size()) 
	{
		std::cerr<<"one size is "<<this->size()<<" and another size is "<<vec.size()<<std::endl;
		throw COException("VectorBase dot operation error: the length of two VectorBases do not equal to each other");
	}
	else{
		scalar sum = blas::copt_blas_dot(this->size(),this->dataPtr(),this->interval(),vec.dataPtr(),vec.interval());
		return sum;
	}
}

template<class scalar,class index>
void VectorBase<scalar,index>::scale(scalar s)
{
	blas::copt_blas_scal(this->size(),s,this->dataPtr(),this->interval());
}

template<class scalar,class index>
VectorBase<scalar,index> VectorBase<scalar,index>::operator* (scalar s)
{
	VectorBase result;
	result.setArray(this->size(),this->dataPtr(),this->interval());
	result.scale(s);
	return result;
}

template<class scalar,class index>
VectorBase<scalar,index> VectorBase<scalar,index>::operator+ (const VectorBase& vec) const
{
	if(this->size()!=vec.size())
	{
		std::cerr<<"one size is "<<this->size()<<" another size is "<<vec.size()<<std::endl;
		throw COException("VectorBase summation error: the length of two VectorBases do not equal to each other");
	}
	VectorBase result(this->size());
	for ( index i = 0 ; i < this->size() ; ++ i ){
		result[i] = this->operator[](i)+vec[i];
	}
	return result;
}

template<class scalar,class index>
VectorBase<scalar,index> VectorBase<scalar,index>::operator- (const VectorBase& vec) const
{
	if(this->size()!=vec.size()) 
	{
		std::cerr<<"one size is "<<this->size()<<" another size is "<<vec.size()<<std::endl;
		throw COException("VectorBase summation error: the length of two VectorBases do not equal to each other");
	}
	VectorBase result(this->size());
	for ( index i = 0 ; i < this->size() ; ++ i ){
		result[i] = this->operator[](i)-vec[i];
	}
	return result;
}

template<class scalar,class index>
VectorBase<scalar,index> VectorBase<scalar,index>::operator- () const
{
	VectorBase result(this->size());
	for ( index i = 0 ; i < this->size() ; ++ i ){
		result[i] = -this->operator[](i);
	}
	return result;
}

template<class scalar,class index>
VectorBase<scalar,index> VectorBase<scalar,index>::vecE(index size, index i )
{
	if ( i < 0 || i >= size ) throw COException("Index error: out of range!");
	VectorBase vec(size);
	vec[i] = 1.0;
	return vec;
}

template<class scalar,class index>
VectorBase<scalar,index> VectorBase<scalar,index>::vecE(index size,index i,const scalar s)
{
	if ( i < 0 || i >= size ) throw COException("Index error: out of range!");
	VectorBase vec(size);
	vec[i] = s;
	return vec;
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
	for ( const auto& s : indices ){
		result[i] = this->operator[](s);
		++ i;
	} 
	return result;
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
template<class InputIterator>
VectorBase<scalar,index> VectorBase<scalar,index>::block(const index bs, InputIterator begin, InputIterator end)const
{
	VectorBase result(bs);
	index i = 0;
	for ( auto s = begin ; s != end ; ++ s )
	{
		result[i] = this->operator[](*s);
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
	for ( const auto& s : indices ){
		this->operator[](i) = vec[s];
		++ i;
	}
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
template<class InputIterator>
void VectorBase<scalar,index>::blockFromVector(const VectorBase& vec, const index bs, InputIterator begin, InputIterator end)
{
	this->resize(bs);
	index i = 0;
	for ( auto iter = begin ; iter != end ; ++ iter )
	{
		this->operator[](i) = vec[*iter];
		++i;
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

template<class scalar,class index>
VectorBase<scalar,index> VectorBase<scalar,index>::random( const index s )
{
	std::uniform_real_distribution<typename get_pod_type<scalar>::type> unif(0.0,1.0);
	VectorBase result(s);
	for ( int i = 0 ; i < s; ++ i )
	{
		if(is_real<scalar>::value)
			ForceAssignment(unif(copt_rand_eng),result(i));
		else
			ForceAssignment(std::complex<typename get_pod_type<scalar>::type>(unif(copt_rand_eng),unif(copt_rand_eng)),result(i));
	}
	return result;
}
}// End of namespace COPT


#endif