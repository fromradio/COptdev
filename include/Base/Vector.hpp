//		This file is part of open library COPT
//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef VECTOR_H
#define VECTOR_H
/*
	basic operation that will be used in the pipeline
*/
/*
	the exception class for COPT
*/

namespace COPT
{
// declaration
template <class FT>
class MatrixBase;

template <class FT>
class VectorBase : public Array<FT>{
public:
	// the type of float number
	typedef 				Array<FT>				Arr;
	typedef 				FT						ScalarType;

public:

	// Constructor

	/*
		default constructor
	*/
	VectorBase();
	/*
		Construct the vector with specific length
			if data is NULL, a zero vector is constructed
	*/

	VectorBase( int length , ScalarType* data = NULL )
		:
		Arr(length,data)
	{
	}

	/*
		Copy assignment
	*/
	VectorBase ( const VectorBase& vec )
		
	{
		this->__size=vec.size();
		this->__data_ptr=new ScalarType[vec.size()];
		blas::copt_blas_copy(this->__size,vec.dataPtr(),1,this->__data_ptr,1);
	}

	/*
		API with vector in stdlib
	*/

	VectorBase(const std::vector<ScalarType>& vec)
	{
		this->__size=static_cast<size_t>(vec.size());
		this->__data_ptr=new ScalarType [vec.size()];
		for ( int i = 0 ; i < this->__size ; ++ i )
			this->__data_ptr[i] = vec[i];
	}

#ifdef EIGEN
	VectorBase(const Eigen::Matrix<ScalarType,Eigen::Dynamic,1>& vec)
		:
		Arr(vec.size())
	{
		for ( int i = 0 ; i < this->__size ; ++ i ){
			this->__data_ptr[i]  =vec(i);
		}
	}
#endif

	/*
		Deconstructor
	*/

	~VectorBase()
	{
	}

	/*
		Copyt operation
	*/

	VectorBase& operator= (const VectorBase& vec ){
		this->resize(vec.size());
		this->setData(vec.size(),vec.dataPtr());
		return *this;
	}

	/** Matlab-like element assignment */
	ScalarType& operator() (const size_t i );
	const ScalarType& operator() (const size_t i )const ;

	/** overload operations*/
	//%{
	/** operator< */
	bool operator< (const VectorBase& vec)const;

	/** operator<= */
	bool operator<=(const VectorBase& vec)const;

	/** operator> */
	bool operator> (const VectorBase& vec)const;

	/** operator>= */
	bool operator>=(const VectorBase& vec)const;

	/** operator== */
	bool operator==(const VectorBase& vec)const;

	/** operator!= */
	bool operator!=(const VectorBase& vec)const;
	//%}

	/*
		Mathematical operations
	*/
	/*
	 * 			Square norm of the VectorBase
	 */
	ScalarType squaredNorm() const{
		ScalarType result = 0;
		for ( int i = 0 ; i < this->__size ; ++ i ){
			result += this->__data_ptr[i]*this->__data_ptr[i];
		}
		return result;
	}
	// dot operation
	ScalarType dot(const VectorBase& vec) const{
		if(this->__size!=vec.size()) throw COException("VectorBase error: the length of two VectorBases do not equal to each other");
		else{
			ScalarType sum = blas::copt_blas_dot(this->__size,this->__data_ptr,1,vec.dataPtr(),1);
			return sum;
		}
	}

	// scale with a special length
	void scale(ScalarType s){
		blas::copt_blas_scal(this->__size,s,this->__data_ptr,1);
	}
	// multiply with a scalar
	VectorBase operator* (ScalarType s){
		VectorBase result(*this);
		result.scale(s);
		return result;
	}

	friend VectorBase operator* (ScalarType s,const VectorBase& vec){
		VectorBase result(vec);
		result.scale(s);
		return result;
	}

	// summation operation
	VectorBase operator+ (const VectorBase& vec) const{
		if(this->__size!=vec.size()) throw COException("VectorBase error: the length of two VectorBases do not equal to each other");
		VectorBase<ScalarType> result(this->__size);
		for ( int i = 0 ; i < this->__size ; ++ i ){
			result[i] = this->__data_ptr[i]+vec[i];
		}
		return result;
	}

	//subtraction operation
	VectorBase operator- (const VectorBase& vec) const{
		if(this->__size!=vec.size()) throw COException("VectorBase error: the length of two VectorBases do not equal to each other");
		VectorBase<ScalarType> result(this->__size);
		for ( int i = 0 ; i < this->__size ; ++ i ){
			result[i] = this->__data_ptr[i]-vec[i];
		}
		return result;
	}

	// 
	VectorBase operator- () const{
		VectorBase<ScalarType> result(this->__size);
		for ( int i = 0 ; i < this->__size ; ++ i ){
			result[i] = -this->__data_ptr[i];
		}
		return result;
	}

	/* overload of stream
	 */
	friend std::ostream& operator<<(std::ostream& os,const VectorBase& vec){
		os<<"[ ";
		for ( int i = 0 ; i< vec.size()-1 ; ++ i ){
			os<<vec[i]<<" , ";
		}
		os<<vec[vec.size()-1]<<" ]";
		return os;
	}


	/*		Generate special VectorBases
	 */


	/*			Generate e VectorBase = [0,0,...,1,0,...]
	 *			/param size:		the size of the VectorBase
	 *			/param i:			the index of non-zero element
	 */
	static VectorBase vecE(size_t size,int i)
	{
		if ( i < 0 || i >= size ) throw COException("Index error: out of range!");
		VectorBase vec(size);
		vec[i] = 1.0;
		return vec;
	}
	/*			Generate e VectorBase = [0,0,...,s,0,...]
	 *			/param size:		the size of the VectorBase
	 *			/param i:			the index of non-zero element
	 *			/param s:			the value of non-zero element
	 */
	static VectorBase vecE(size_t size,int i,const ScalarType s)
	{
		if ( i < 0 || i >= size ) throw COException("Index error: out of range!");
		VectorBase vec(size);
		vec[i] = s;
		return vec;
	}

	/*			transpose operations
	 */
	MatrixBase<ScalarType> mulTrans(const VectorBase& vec) const;

	/*			the transpose of the vector multiplies a matrix
	 */
	VectorBase transMul(const MatrixBase<ScalarType>& mat) const;

	/** blocking operations */
	//%{
	VectorBase block(const std::set<size_t>& indices)const;
	void blockFromVector(const VectorBase& vec,const std::set<size_t>& indices);
	VectorBase block(const std::vector<size_t>& indices) const;
	void blockFromVector(const VectorBase& vec,const std::vector<size_t>& indices);
	//%}

	/** combination operations */
	//%{
	/** combination of two vectors */
	void combine(
		const VectorBase& v1,
		const VectorBase& v2);
	/** combination of two vectors taking output as parameter */
	static inline void stCombine(
		const VectorBase& v1,
		const VectorBase& v2,
		VectorBase& v);
	//%}

};

}	// end of namespace COPT

#endif