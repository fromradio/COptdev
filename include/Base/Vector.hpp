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
	/** 	size type 		*/
	typedef 				I 						index;
	/**		define the category 	*/
	typedef 				vector_tag 				Category;
	/**		define the kernel 		*/
	typedef 				KernelTrait<FT,index>	Kernel;

private:
	/**		definitions used in implementation */
	typedef 				Array<FT,index>			Arr;
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

	VectorBase( const index size , scalar* data = NULL )
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
			this->setReferredArray(vec.size(),const_cast<scalar*>(vec.dataPtr()),vec.interval());
		}
		else{
			// copy the vector data
			this->setArray(vec.size(),vec.dataPtr(),vec.interval());
		}
	}

	VectorBase( const index size , const referred_array& tag , scalar* data ,const index inter = 1)
		:
		Arr(size,tag,data,inter)
	{
	}

	VectorBase( const index size , const referred_array& tag , const scalar* data ,const index inter = 1)
		:
		Arr(size,tag,const_cast<scalar*>(data),inter)
	{
	}

	/*
		API with vector in stdlib
	*/

	VectorBase(const std::vector<scalar>& vec)
		:
		Arr()
	{
		this->resize(vec.size(),1);
		for ( index i = 0 ; i < this->size() ; ++ i )
			this->operator[](i) = vec[i];
	}

#ifdef EIGEN
	VectorBase(const Eigen::Matrix<scalar,Eigen::Dynamic,1>& vec)
		:
		Arr(vec.size())
	{
		for ( index i = 0 ; i < this->size() ; ++ i ){
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
			this->setReferredArray(vec.size(),const_cast<scalar*>(vec.dataPtr()),vec.interval());
		}
		else{
			this->setArray(vec.size(),vec.dataPtr(),vec.interval());
		}
		return *this;
	}

	/** Matlab-like element assignment */
	scalar& operator() (const index i );
	const scalar& operator() (const index i )const ;

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
	scalar squaredNorm() const{
		scalar result = 0;
		for ( index i = 0 ; i < this->size() ; ++ i ){
			result += this->operator[](i)*this->operator[](i);
		}
		return result;
	}
	/** normalize current vector and previous norm is returned*/
	scalar normalize();
	
	// dot operation
	scalar dot(const VectorBase& vec) const{
		if(this->size()!=vec.size()) throw COException("VectorBase error: the length of two VectorBases do not equal to each other");
		else{
			scalar sum = blas::copt_blas_dot(this->size(),this->dataPtr(),this->interval(),vec.dataPtr(),vec.interval());
			return sum;
		}
	}

	// scale with a special length
	void scale(scalar s){
		blas::copt_blas_scal(this->size(),s,this->dataPtr(),this->interval());
	}
	// multiply with a scalar
	VectorBase operator* (scalar s){
		VectorBase result;
		result.setArray(this->size(),this->dataPtr(),this->interval());
		result.scale(s);
		return result;
	}

	friend VectorBase operator* (scalar s,const VectorBase& vec){
		VectorBase result;
		result.setArray(vec.size(),const_cast<scalar*>(vec.dataPtr()),vec.interval());
		result.scale(s);
		return result;
	}

	// summation operation
	VectorBase operator+ (const VectorBase& vec) const{
		if(this->size()!=vec.size()) throw COException("VectorBase error: the length of two VectorBases do not equal to each other");
		VectorBase result(this->size());
		for ( index i = 0 ; i < this->size() ; ++ i ){
			result[i] = this->operator[](i)+vec[i];
		}
		return result;
	}

	//subtraction operation
	VectorBase operator- (const VectorBase& vec) const{
		if(this->size()!=vec.size()) throw COException("VectorBase error: the length of two VectorBases do not equal to each other");
		VectorBase result(this->size());
		for ( index i = 0 ; i < this->size() ; ++ i ){
			result[i] = this->operator[](i)-vec[i];
		}
		return result;
	}

	// 
	VectorBase operator- () const{
		VectorBase result(this->size());
		for ( index i = 0 ; i < this->size() ; ++ i ){
			result[i] = -this->operator[](i);
		}
		return result;
	}

	/* overload of stream
	 */
	friend std::ostream& operator<<(std::ostream& os,const VectorBase& vec){
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


	/*		Generate special VectorBases
	 */


	/*			Generate e VectorBase = [0,0,...,1,0,...]
	 *			/param size:		the size of the VectorBase
	 *			/param i:			the index of non-zero element
	 */
	static VectorBase vecE(index size,index i)
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
	static VectorBase vecE(index size,index i,const scalar s)
	{
		if ( i < 0 || i >= size ) throw COException("Index error: out of range!");
		VectorBase vec(size);
		vec[i] = s;
		return vec;
	}

	/*			transpose operations
	 */
	MatrixBase<scalar,index> mulTrans(const VectorBase& vec) const;

	/*			the transpose of the vector multiplies a matrix
	 */
	VectorBase transMul(const MatrixBase<scalar,index>& mat) const;

	/** blocking operations */
	//%{
	VectorBase block(const std::set<index>& indices)const;
	void blockFromVector(const VectorBase& vec,const std::set<index>& indices);
	VectorBase block(const std::vector<index>& indices) const;
	void blockFromVector(const VectorBase& vec,const std::vector<index>& indices);
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