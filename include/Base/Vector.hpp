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


#ifndef VECTOR_HPP__
#define VECTOR_HPP__
/*
	basic operation that will be used in the pipeline
*/
/*
	the exception class for COPT
*/

namespace COPT
{
// declaration
template <class FT,class I>
class MatrixBase;

template <class FT,class I = int>
class VectorBase 
	: 
	public Array<FT,I>
{
public:
	/** 	scalar type 	*/
	typedef 				FT						scalar;
	/** 	the pod type of scalar */
	typedef typename get_pod_type<scalar>::type 	podscalar;
	/** 	size type 		*/
	typedef 				I 						index;
	/**		define the category 	*/
	typedef 				vector_object 			ObjectCategory;
	/**		define the kernel 		*/
	typedef 				KernelTrait<FT,index>	Kernel;

private:
	/**		definitions used in implementation */
	typedef 				Array<FT,index>			Arr;
public:

	
	/** constructors and deconstructor */
	//%{
	/** default constructor */
	VectorBase();
	
	/*
		Construct the vector with specific length
			if data is NULL, a zero vector is constructed
	*/
	VectorBase(const index size, scalar* data = nullptr);

	/** Copy assignment */
	VectorBase (const VectorBase& vec);

	VectorBase(const index size, const referred_array& tag, scalar* data, const index inter = 1);

	VectorBase(const index size ,const referred_array& tag, const scalar* data, const index inter = 1);

	/** API with vector in standard library */
	VectorBase(const std::vector<scalar>& vec);

	/** Deconstructor */
	~VectorBase();
	//%}

	/** copy operation */
	VectorBase& operator=(const VectorBase& vec );

	/** Matlab-like element assignment */
	scalar& operator()(const index i );
	const scalar& operator()(const index i )const ;

	/** overload operations*/
	//%{
	/** operator< */
	bool operator< (const VectorBase& vec)const;
	bool operator< (const scalar s)const;

	/** operator<= */
	bool operator<=(const VectorBase& vec)const;
	bool operator<=(const scalar s)const;

	/** operator> */
	bool operator> (const VectorBase& vec)const;
	bool operator> (const scalar s)const;

	/** operator>= */
	bool operator>=(const VectorBase& vec)const;
	bool operator>=(const scalar s)const;

	/** operator== */
	bool operator==(const VectorBase& vec)const;
	bool operator==(const scalar s)const;

	/** operator!= */
	bool operator!=(const VectorBase& vec)const;
	bool operator!=(const scalar s)const;
	//%}

	/*
		Mathematical operations
	*/
	/*
	 * 			Square norm of the VectorBase
	 */
	podscalar squaredNorm() const;

	/** l1 norm */
	podscalar absNorm() const;

	/** l2 norm */
	auto norm() const->podscalar;

	/** normalize current vector and previous norm is returned*/
	podscalar normalize();
	
	/** dot operation */
	scalar dot(const VectorBase& vec) const;

	/** scale with a special length */
	void scale(scalar s);

	/** multiply with a scalar */
	VectorBase operator* (scalar s);


	friend VectorBase operator* (scalar s, const VectorBase& vec){
		VectorBase result;
		result.setArray(vec.size(),const_cast<scalar*>(vec.dataPtr()),vec.interval());
		result.scale(s);
		return result;
	}

	/** summation operation */
	VectorBase operator+ (const VectorBase& vec) const;

	/** subtraction operation */
	VectorBase operator- (const VectorBase& vec) const;

	/** */
	VectorBase operator- () const;

	/** overload of output stream */
	friend std::ostream& operator<<(std::ostream& os, const VectorBase& vec){
		if ( vec.size() == 0 )
		{
			os<<"[ ]";
			return os;
		}
		os<<"[ ";
		for ( index i = 0 ; i< vec.size()-1 ; ++ i ){
			os<<vec[i]<<" , ";
		}
		os<<vec[vec.size()-1];
		os<<" ]";
		return os;
	}

	/** Generate special VectorBases */
	//%{
	/*			Generate e VectorBase = [0,0,...,1,0,...]
	 *			/param size:		the size of the VectorBase
	 *			/param i:			the index of non-zero element
	 */
	static VectorBase vecE(index size, index i);

	/*			Generate e VectorBase = [0,0,...,s,0,...]
	 *			/param size:		the size of the VectorBase
	 *			/param i:			the index of non-zero element
	 *			/param s:			the value of non-zero element
	 */
	static VectorBase vecE(index size, index i, const scalar s);
	//%}

	/** transpose operations */
	MatrixBase<scalar,index> mulTrans(const VectorBase& vec) const;

	/*			the transpose of the vector multiplies a matrix
	 */
	VectorBase transMul(const MatrixBase<scalar,index>& mat) const;

	/** blocking operations */
	//%{
	VectorBase block(const std::set<index>& indices)const;
	VectorBase block(const std::vector<index>& indices) const;
	template<class InputIterator>
	VectorBase block(const index bs, InputIterator begin, InputIterator end)const;
	void blockFromVector(const VectorBase& vec, const std::set<index>& indices);
	void blockFromVector(const VectorBase& vec, const std::vector<index>& indices);
	template<class InputIterator>
	void blockFromVector(const VectorBase& vec, const index bs, InputIterator begin, InputIterator end);
	//%}

	/** combination operations */
	//%{
	/** combination of two vectors */
	void combine(const VectorBase& v1, const VectorBase& v2);
	/** combination of two vectors taking output as parameter */
	static inline void stCombine(const VectorBase& v1, const VectorBase& v2, VectorBase& v);
	//%}

	/** generate a random vector */
	static inline VectorBase random(const index i);

};

}	// end of namespace COPT

#endif