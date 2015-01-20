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


#ifndef MATRIX_BASE_HPP__
#define MATRIX_BASE_HPP__

/*
 *				class 'MatrixBase'
 *					dervied from base class 'Array'
 */
namespace COPT
{
/*
 *	Class of 'MatrixBase'
 *		the data is stored column by column
 *		the matrix is assumed not to be transpose or symmetric
 *		but symmetric and transpose flag is designed
 */
template<class FT,class I = int,int RowAtCompileTime=Dynamic,int ColAtCompileTime=Dynamic>
class MatrixBase
	:
	public Array<FT,I,(RowAtCompileTime==Dynamic||ColAtCompileTime==Dynamic)?Dynamic:(RowAtCompileTime+1)*(ColAtCompileTime+1)>
{
public:

	/** the scalar type */
	typedef 			FT 				 			scalar;
	/** pod scalar type */
	typedef typename get_pod_type<scalar>::type 	podscalar;
	/** the size type used */
	typedef 			I 							index;						
	/** define fthe category */
	typedef 			matrix_object				ObjectCategory;
	/** define the trait */
	typedef 			KernelTrait<FT,index>		Kernel;

private:
	
	/** definition used in implementation */
	typedef 			VectorBase<FT,index>			Vector;
	typedef 			COPT::Array<FT,index>			Array;

	/**			private variables			*/
	//%{
	/** the size of rows */
	index					__rows;
	/** the size of columns */
	index 					__cols;
	/** whether the matrix is symmetric */
	bool					__sym;
	/** whether the matrix is transpose */
	bool 					__trans;
	/** lda */
	index 					__lda;
	/** is row number dynamic */
	bool 					__is_row_dynamic;
	/** is column number dynamic */
	bool 					__is_col_dynamic;
	//%}

public:
	/** constructor and deconstructor */
	//%{
	/** default constructor */
	MatrixBase();

	MatrixBase(const index m, const index n, const scalar *data=NULL);

	/** Copy assignment */
	MatrixBase(const MatrixBase& mat);

	/** deconstructor */
	~MatrixBase();
	//%} end of constructor and deconstructor



	/**	getters and setters*/
	//%{
	/** get the number of rows */
	const index&		rows() const;
	/** is the row number dynamic */
	bool isRowDynamic() const;
	/** get the number of columns */
	const index& 		cols() const;
	bool isColumnDynamic() const;

	/**	matlab-like element getter */
	scalar& operator() (const index i, const index j);

	const scalar& operator() (const index i, const index j) const;

	/** get the element using Arr */
	const scalar& data(const index i) const;

	/**	obtain the i-th row */
	Vector row(const index num);
	const Vector row(const index num) const;
	/**	obtain the i-th column */
	Vector col(const index num);
	const Vector col(const index num) const;


	/** set element using Arr */

	void set (const index i, const scalar value);


	/** resize the matrix */
	void resize (index m, index n);

	/** set the matrix to be symmetric */
	void setSymmetricFlag(bool sym);
	bool isSymmetric() const;

	/** set the matrix to be transpose mode */
	void setTransposeFlag(bool sym);
	bool isTranspose() const;

	/** Copy operation */
	MatrixBase& operator= (const MatrixBase& mat);

	/*
		Mathematical operations
	*/

	/** summation */
	MatrixBase operator+ (const MatrixBase& mat);

	// subtraction
	// need to be tested
	MatrixBase operator- (const MatrixBase& mat);

	/** matrix multiplications */
	VectorBase<scalar,index> operator* (const VectorBase<scalar,index>& vec) const;

	/** matrix and matrix multiplication */
	MatrixBase operator* (const MatrixBase& mat) const;


	// multiplication between a scalar and a matrix
	friend MatrixBase operator* (const scalar s, const MatrixBase& mat)
	{
		MatrixBase result(mat.rows(),mat.cols());
		for ( index i = 0 ; i < mat.rows() ; ++ i )
			for ( index j = 0 ; j < mat.cols() ; ++ j )
				result(i,j) = mat(i,j)*s;
		return result;
	}

	// multiplication between a scalar and a matrix
	friend MatrixBase operator* (const MatrixBase& mat, const scalar s)
	{
		return s*mat;
	}

	// transpose
	MatrixBase transpose() const;

	/** transpose muliplication */
	Vector transMulti(const Vector& vec) const;

	MatrixBase transMulti(const MatrixBase& mat) const;

	/*				overload of ostream
	 */
	friend std::ostream& operator<<(std::ostream& os, const MatrixBase& mat)
	{
		for ( index i = 0 ; i < mat.rows() ; ++ i ){
			for ( index j = 0 ; j < mat.cols() ; ++ j )
				os<<mat(i,j)<<' ';
			os<<std::endl;
		}
		return os;
	}

	/*				solve linear system
	 *
	 */
	 Vector lapackSolve(const Vector& vec){
		if (__rows != __cols )
			throw COException("Solving Error: the matrix is not square!");
		if (__rows != vec.size() )
			throw COException("Solving Error: the size of right hand vector is wrong!");
		Vector v(vec);
		int pivot[__rows], c2 = 1,info;
		scalar* a = new scalar[__rows*__cols];
		blas::copt_blas_copy(__rows*__cols,this->dataPtr(),1,a,1);
		copt_lapack_gesv(&__rows,&c2,a,
			&__rows,pivot,v.dataPtr(),
			&__rows,
			&info);
		delete[] a; // delete the temporaray array
		return v;
	}
	
	Vector solve(const Vector& vec){
		return lapackSolve(vec);
	}

	Vector leastSquareSolve(const Vector& vec){
		if( __cols != vec.size() )
			throw COException("least square error: the size is not consistent!");
		scalar *a = new scalar[__rows*__cols];
		scalar *b = new scalar[vec.size()];
		scalar *s = new scalar[std::min(__rows,__cols)];
		index info = -1;
		index rank;
		blas::copt_blas_copy(__rows*__cols,this->dataPtr(),1,a,1);
		blas::copt_blas_copy(vec.size(),vec.dataPtr(),1,b,1);
		copt_lapack_gelss(__rows,__cols,1,a,__rows,b,vec.size(),s,-1.0,&rank,&info);
		Vector result(__cols,b);
		delete[]a;
		delete[]b;
		return result;
	}
	


	/*			Special MatrixBase
	 *
	 */
	static MatrixBase identity(index m, index n){
		MatrixBase result(m,n);
		// std::cout<<result<<std::endl;
		index min = std::min(m,n);
		for ( index i = 0 ; i < min ; ++ i )
			result(i,i) = static_cast<scalar>(1.0);
		return result;
	}

	static MatrixBase identity(index m, index n, const scalar s);

	/** Blocking methods */
	//%{
	
	/** Blocking matrix from a given matrix with specific row numbers and column numbers*/
	void blockFromMatrix(
		const MatrixBase& mat,
		const std::set<index>& rownums,
		const std::set<index>& colnums);
	
	/** Blocking matrix from a given matrix with just columns */
	void columnBlockFromMatrix(
		const MatrixBase& mat,
		const std::set<index>& colnums);
	
	/** Blocking matrix from a given matrix with just rows */
	void rowBlockFromMatrix(
		const MatrixBase& mat,
		const std::set<index>& rownums);

	/** Blocking matrix from a given matrix, order is not considered */
	void blockFromMatrix(
		const MatrixBase& mat,
		const std::vector<index>& rownums,
		const std::vector<index>& colnums);

	/** blocking matrix using iterator */
	template<class InputIterator>
	void blockFromMatrix(
		const MatrixBase& mat,
		const InputIterator& rowbegin,
		const InputIterator& rowend,
		const InputIterator& colbegin,
		const InputIterator& colend );

	/** Blocking matrix from a given matrix, order is not considered */
	void columnBlockFromMatrix(
		const MatrixBase& mat,
		const std::vector<index>& colnums);

	/** column blocking from a given matrix */
	template<class InputIterator>
	void columnBlockFromMatrix(
		const MatrixBase& mat,
		const InputIterator& colbegin,
		const InputIterator& colend);

	/** Blocking matrix from a given matrix, order is not considered */
	void rowBlockFromMatrix(
		const MatrixBase& mat,
		const std::vector<index>& rownums);

	template<class InputIterator>
	void rowBlockFromMatrix(
		const MatrixBase& mat,
		const InputIterator& rowbegin,
		const InputIterator& rowend);

	//%}

	/** combining methods */
	//%{

	/** combining along row direction */
	void combineAlongRow(
		const MatrixBase& m1,
		const MatrixBase& m2);
	/** combining along column direction */
	void combineAlongColumn(
		const MatrixBase& m1,
		const MatrixBase& m2);

	/** combine along row direction taking matrix as parameter */
	static inline void stCombineAlongRow(
		const MatrixBase& m1,
		const MatrixBase& m2,
		MatrixBase& m);
	/** combine along column direction taking matrix as parameter */
	static inline void stCombineAlongColumn(
		const MatrixBase&m1,
		const MatrixBase& m2,
		MatrixBase& m);
	//%}

	/** set a random matrix */
	void setRandom(const index rows, const index cols);
	static inline MatrixBase random(const index rows, const index cols);
	/** compute A^TA of a given matrix */
	void mtm(MatrixBase& mat) const;

	/** norms */
	//%{
	/** compute the operation norm */
	podscalar operationNorm() const;
	//%}

};// End of class MatrixBase




}// End of namespace COPT
#endif