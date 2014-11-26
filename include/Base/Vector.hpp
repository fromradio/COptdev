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
template <class FT,class Size>
class MatrixBase;

template <class FT,class Size = size_t>
class VectorBase : public Array<FT,Size>{
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

	VectorBase( const Size size , ScalarType* data = NULL )
		:
		Arr(size , data)
	{
	}

	/** Copy assignment */
	VectorBase ( const VectorBase& vec )
		:
		Arr()
	{
		if (vec.isReferred())
		{
			// the vector is a referred vector
			this->setReferredArray(vec.size(),const_cast<ScalarType*>(vec.dataPtr()),vec.interval());
		}
		else{
			// copy the vector data
			this->setArray(vec.size(),vec.dataPtr(),vec.interval());
		}
	}

	VectorBase( const Size size , const referred_array& tag , ScalarType* data ,const Size inter = 1)
		:
		Arr(size,tag,data,inter)
	{
	}

	/*
		API with vector in stdlib
	*/

	VectorBase(const std::vector<ScalarType>& vec)
		:
		Arr()
	{
		this->resize(vec.size(),1);
		for ( int i = 0 ; i < this->size() ; ++ i )
			this->operator[](i) = vec[i];
	}

#ifdef EIGEN
	VectorBase(const Eigen::Matrix<ScalarType,Eigen::Dynamic,1>& vec)
		:
		Arr(vec.size())
	{
		for ( int i = 0 ; i < this->size() ; ++ i ){
			this->operator[](i)  = vec(i);
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
		if (vec.isReferred()){
			this->setReferredArray(vec.size(),const_cast<ScalarType*>(vec.dataPtr()),vec.interval());
		}
		else{
			this->setArray(vec.size(),vec.dataPtr(),vec.interval());
		}
		return *this;
	}

	/** Matlab-like element assignment */
	ScalarType& operator() (const Size i );
	const ScalarType& operator() (const Size i )const ;

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
		for ( int i = 0 ; i < this->size() ; ++ i ){
			result += this->operator[](i)*this->operator[](i);
		}
		return result;
	}
	// dot operation
	ScalarType dot(const VectorBase& vec) const{
		if(this->size()!=vec.size()) throw COException("VectorBase error: the length of two VectorBases do not equal to each other");
		else{
			ScalarType sum = blas::copt_blas_dot(this->size(),this->dataPtr(),this->interval(),vec.dataPtr(),vec.interval());
			return sum;
		}
	}

	// scale with a special length
	void scale(ScalarType s){
		blas::copt_blas_scal(this->size(),s,this->dataPtr(),this->interval());
	}
	// multiply with a scalar
	VectorBase operator* (ScalarType s){
		VectorBase result;
		result.setArray(this->size(),this->dataPtr(),this->interval());
		result.scale(s);
		return result;
	}

	friend VectorBase operator* (ScalarType s,const VectorBase& vec){
		VectorBase result;
		result.setArray(vec.size(),const_cast<ScalarType*>(vec.dataPtr()),vec.interval());
		result.scale(s);
		return result;
	}

	// summation operation
	VectorBase operator+ (const VectorBase& vec) const{
		if(this->size()!=vec.size()) throw COException("VectorBase error: the length of two VectorBases do not equal to each other");
		VectorBase<ScalarType> result(this->size());
		for ( int i = 0 ; i < this->size() ; ++ i ){
			result[i] = this->operator[](i)+vec[i];
		}
		return result;
	}

	//subtraction operation
	VectorBase operator- (const VectorBase& vec) const{
		if(this->size()!=vec.size()) throw COException("VectorBase error: the length of two VectorBases do not equal to each other");
		VectorBase<ScalarType> result(this->size());
		for ( int i = 0 ; i < this->size() ; ++ i ){
			result[i] = this->operator[](i)-vec[i];
		}
		return result;
	}

	// 
	VectorBase operator- () const{
		VectorBase<ScalarType> result(this->size());
		for ( int i = 0 ; i < this->size() ; ++ i ){
			result[i] = -this->operator[](i);
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
	static VectorBase vecE(Size size,int i)
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
	static VectorBase vecE(Size size,int i,const ScalarType s)
	{
		if ( i < 0 || i >= size ) throw COException("Index error: out of range!");
		VectorBase vec(size);
		vec[i] = s;
		return vec;
	}

	/*			transpose operations
	 */
	MatrixBase<ScalarType,Size> mulTrans(const VectorBase& vec) const;

	/*			the transpose of the vector multiplies a matrix
	 */
	VectorBase transMul(const MatrixBase<ScalarType,Size>& mat) const;

	/** blocking operations */
	//%{
	VectorBase block(const std::set<Size>& indices)const;
	void blockFromVector(const VectorBase& vec,const std::set<Size>& indices);
	VectorBase block(const std::vector<Size>& indices) const;
	void blockFromVector(const VectorBase& vec,const std::vector<Size>& indices);
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