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


#ifndef EIGEN_WRAPPER_HPP__
#define EIGEN_WRAPPER_HPP__


namespace COPT
{

/*			A wrapper from open source library 'Eigen's linear solver.
 *			The ldlt solver is chosen since it is the fastest one due to
 *			the description of 'Eigen'.
 */
template<class Matrix>
class EigenSolver
	:
	public LinearSolver<Matrix>
{
private:
	typedef typename Matrix::scalar 				scalar;
	typedef typename Matrix::index 					index;
	typedef VectorBase<scalar,index>				Vector;


	void doCompute(const Matrix& mat){}
	Vector doSolve( const Vector& b ){return Vector();}
	Matrix doSolve( const Matrix& b ){return Matrix();}

public:
	
	EigenSolver(){}
	~EigenSolver(){}

	void clear(){}
};


}

#endif