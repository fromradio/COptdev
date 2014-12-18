//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef CHOLESKY_SOLVER_HPP__
#define CHOLESKY_SOLVER_HPP__

namespace COPT
{

template<class Matrix>
class CholeskySolver
{
private:
	typedef typename Matrix::index 				index;
	typedef typename Matrix::scalar 			scalar;
	typedef VectorBase<scalar,index>			Vector;

	/** the data */
	scalar 							*__a;

	/** the order of matrix */
	index 							__n;

	/** lda */
	index 							__lda;

	/** info */
	index 							__info;

	/** since the input matrix has to be symmetric semidefinite */
	bool validation( const Matrix& mat ) const;

public:

	/** constructor and deconstructor */
	//%{
	/** default constructor */
	CholeskySolver();
	CholeskySolver(const Matrix& mat);
	~CholeskySolver();
	//%}

	/** compute the problem */
	void compute( const Matrix& mat );

	/** solving the problem */
	Vector solve ( const Vector& b );

	/** solving the problem */
	Matrix solve ( const Matrix& b );

	/** inver matrix */
	Matrix inverse ();

	/** validation */
	void squareValidation( ) const;

	/** clear */
	void clear();
};

template<class Matrix>
CholeskySolver<Matrix>::CholeskySolver()
	:
	__a(NULL)
{
}

template<class Matrix>
CholeskySolver<Matrix>::CholeskySolver( const Matrix& mat )
	:
	__a(NULL)
{
	compute(mat);
}

template<class Matrix>
CholeskySolver<Matrix>::~CholeskySolver()
{
	clear();
}

template<class Matrix>
void CholeskySolver<Matrix>::compute( const Matrix& mat )
{
	if (!validation(mat))
	{
		std::cerr<<"Warning in cholesky solver: the input matrix is not symmetric. Please confirm that!"<<std::endl;
		return;
	}
	SAFE_DELETE_ARRAY(__a);
	__a = new scalar[mat.size()];
	__n = mat.rows();
	__lda = mat.rows();
	blas::copt_blas_copy(mat.size(),mat.dataPtr(),1,__a,1);
	copt_lapack_potrf('U',__n,__a,__lda,&__info);
	if( __info != 0 )
		std::cerr<<"Warning in Cholesky solver: something computation is not wrong!"<<std::endl;
}

template<class Matrix>
typename CholeskySolver<Matrix>::Vector CholeskySolver<Matrix>::solve( const Vector& b )
{
	if ( __n != b.size() )
	{
		std::cerr<<"The order of matrix is "<<__n<<" and the dimension of vector is "<<b.size()<<std::endl;
		throw COException("CholeskySolver solving error: the size if not consistent!");
	}
	Vector result(b);
	copt_lapack_potrs('U',__n,1,__a,__lda,result.dataPtr(),result.size(),&__info);
	if ( __info != 0 )
		std::cerr<<"Warning in Cholesky solver: solving is wrong!"<<std::endl;
	return result;
}

template<class Matrix>
Matrix CholeskySolver<Matrix>::solve( const Matrix& b )
{
	if ( __n != b.rows() )
	{
		std::cerr<<"The order of matrix is "<<__n<<" and the dimension of vector is "<<b.rows()<<std::endl;
		throw COException("CholeskySolver solving error: the size if not consistent!");
	}
	Matrix result(b);
	copt_lapack_potrs('U',__n,b.cols(),__a,__lda,result.dataPtr(),result.rows(),&__info);
	if ( __info != 0 )
		std::cerr<<"Warning in Cholesky solver: solving is wrong!"<<std::endl;
	return result;
}

template<class Matrix>
bool CholeskySolver<Matrix>::validation( const Matrix& mat )const
{
	return mat.isSymmetric();
}

template<class Matrix>
Matrix CholeskySolver<Matrix>::inverse( )
{
	Matrix result(__n,__n,__a);
	copt_lapack_potri('U',__n,result.dataPtr(),__lda,&__info);
	result.setSymmetricFlag(true);
	return result;
}

template<class Matrix>
void CholeskySolver<Matrix>::squareValidation() const
{
}

template<class Matrix>
void CholeskySolver<Matrix>::clear()
{
	SAFE_DELETE_ARRAY(__a);
}

}

#endif