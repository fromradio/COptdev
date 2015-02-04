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


#ifndef ARITH_OPERATION_HPP__
#define ARITH_OPERATION_HPP__

namespace COPT
{

/** check the dimension of two vectors */
template<class Vector>
void checkDimension(const Vector &v1, const Vector &v2)
{
	if(v1.size()!=v2.size())
		throw COException("COPT error: The size of two vectors are not consistent!");
}

/** compute the norm */
template<class T>
typename T::podscalar norm(const T& t,int l)
{
	return norm(t,l,typename T::ObjectCategory());
}

/** compute norm of a vector */
template<class Vector>
typename Vector::podscalar norm(const Vector& vec, int l, const vector_object&)
{
	typedef typename Vector::index ind;
	typedef typename Vector::podscalar podscalar;
	switch(l)
	{
		case 0:
		{
			ind s = vec.size();
			podscalar norm=0;
			for ( ind i = 0 ; i<s ; ++i )
			{
				norm += (IS_ZERO(vec[i]))?0:1;
			}
			return norm;
		}
		break;
		case 1:
		{
			return vec.absNorm();
		}
		break;
		case 2:
		{
			return vec.norm();
		}
		break;
		default:
		{
			throw COException("Unsupported type of norm operation");
			return 0;
		}
		break;
	}
}

/** compute a general norm of a vector */
template<class Vector>
typename Vector::podscalar norm(const Vector &vec, typename Vector::podscalar l, const vector_object&)
{
	
}

/** compute the distance between two vectors */
template<class Vector>
typename Vector::podscalar distance(const Vector &v1, const Vector &v2)
{
	checkDimension(v1,v2);
	return (v1-v2).norm();
}

/** compute the mean vector of a set of given vector */
template<class VectorIterator>
typename VectorIterator::value_type mean(VectorIterator begin,VectorIterator end)
{
	typename VectorIterator::value_type vec((*begin).dimension());
	int n=0;
	std::for_each(begin,end,[&vec,&n](typename VectorIterator::value_type& v){vec=vec+v;++n;});
	vec.scale(1.0/n);
	return vec;
}


template<class T>
T sgn(const T& t)
{
	typedef typename T::scalar scalar;
	T result(t);
	for_each(t.begin(),t.end(),[](scalar& s){s=(s<0)?-1:(s>0?1:0);});
}

/** add sparse noise */
template<class T,class S,class I>
void addSparseNoise(T& t,const I sp,const S n)
{
	std::vector<int> tt(t.size());
	for ( int i = 0 ; i < tt.size() ; ++i ) tt[i]=i;
	std::random_shuffle(tt.begin(),tt.end());
	std::uniform_real_distribution<typename T::podscalar> unif(-1.0,1.0);
	for ( int i = 0 ; i < sp ; ++ i ) t[tt[i]]+=n*unif(copt_rand_eng);
}


}

#endif