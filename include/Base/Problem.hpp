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


#ifndef PROBLEM_HPP__
#define PROBLEM_HPP__

namespace COPT
{

/*		This file introduces the basic problem type in COPT. 
 *
 *
 */
template <class kernel>
class VectorProblem
{

	typedef typename kernel::scalar 		scalar;
	typedef typename kernel::podscalar 		podscalar;
	typedef typename kernel::index 			index;
	typedef typename kernel::Matrix 		Matrix;
	typedef typename kernel::Vector 		Vector;

	/** the dimension of the problem */
	index 					__dim;

public:

	typedef kernel 							KernelTrait;

	/** constructor and deconstructor */
	//%{
	/** default constructor */
	VectorProblem( const index dim = 0 );
	/** deconstructor */
	virtual ~VectorProblem(){}
	//%}

	/** compute the objective function */
	virtual podscalar objective( const Vector& x ) const = 0;

	/** check whether the input is valid */
	virtual bool isValidInput( const Vector& x ) const = 0;

	/** check whether the problem is a valid problem */
	virtual bool isValid( ) const = 0;

	/** validation of the problem */
	virtual void validation() const;

	/** getter and setter */
	//%{
	void setDimension( const index dim );
	index dimension() const;
	//%}
};

/******************************Implementation of 'VectorProblem'**********************/
template<class kernel>
VectorProblem<kernel>::VectorProblem( const index dim )
	:
	__dim(dim)
{
}

template<class kernel>
void VectorProblem<kernel>::validation() const
{
	if(!this->isValid())
	{
		throw COException("Vector problem error: Validation of this problem does not pass!");
	}
}

template<class kernel>
void VectorProblem<kernel>::setDimension( const index dim )
{
	__dim = dim;
}

template<class kernel>
typename VectorProblem<kernel>::index VectorProblem<kernel>::dimension() const
{
	return __dim;
}

}// End of namespace COPT

#endif