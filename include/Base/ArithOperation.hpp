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

}

#endif