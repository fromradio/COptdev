//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU


#ifndef EIGEN_VALUE_HPP__
#define EIGEN_VALUE_HPP__

namespace COPT
{
template<class Matrix>
class SymmetricToTridiagonal
{
private:
	typedef typename Matrix::Kernel 				kernel;
	typedef typename kernel::scalar 				scalar;
	typedef typename kernel::index 					index;

	/** uplo */
	char 			__uplo;
	/** the dimension of symmetric matrix */
	index 			__n;
	/** the diagonal elements */
	scalar 			*__d;
	/** the sub-diagonal elements */
	scalar			*__e;
	/** tau */
	scalar 			*__tau;
	/** information */
	index 			__info;
public:
	/** constructor and deconstructor */
	//%{
	SymmetricToTridiagonal();
	SymmetricToTridiagonal( const Matrix& mat );
	~SymmetricToTridiagonal();
	//%}

	/** kernel function */
	void compute(const Matrix& mat);

	/** getter */
	//%{
	index n() const;
	const scalar* d() const;
	const scalar* e() const;
	//%}
};

template<class Matrix>
class PartialEigenSolver
{
private:
	typedef typename Matrix::Kernel 			kernel;
	typedef typename kernel::index 				index;
	typedef typename kernel::scalar 			scalar;
	typedef SymmetricToTridiagonal<Matrix>		tri_generator;

	/** the triagonal */
	tri_generator 		__tri;
	/** the dimension of the matrix */
	index 				__dim;
	/** the actual number of eigenvalues fount */
	index 				__m;
	/** storing the eigenvalues */
	scalar 				*__w;
	/** block infomation */
	index 				*__iblock;
	/**	nsplit */
	index 				__nsplit;
	/** the split array */
	index 				*__isplit;
	/** the information */
	index 				__info;
	/** whether is solved */
	bool 				__is_solved;

public:
	PartialEigenSolver();
	PartialEigenSolver( const Matrix& mat );
	~PartialEigenSolver();

	void compute(const Matrix& mat);

	void solveIndex( index il, index iu );

	scalar computeLargestEigenvalue();
};

/********************Implementation of 'SymmetricToTridiagonal'****************/

template<class Matrix>
SymmetricToTridiagonal<Matrix>::SymmetricToTridiagonal()
	:
	__uplo('U'),
	__d(NULL),
	__e(NULL),
	__tau(NULL)
{
}

template<class Matrix>
SymmetricToTridiagonal<Matrix>::SymmetricToTridiagonal( const Matrix& mat )
	:
	__uplo('U'),
	__d(NULL),
	__e(NULL),
	__tau(NULL)
{
	compute(mat);
}

template<class Matrix>
void SymmetricToTridiagonal<Matrix>::compute(const Matrix& mat)
{
	if (!mat.isSymmetric())
		std::cerr<<" SymmetricToTridiagonal Warning: Please make sure that whether the input matrix is symmetric or not!"<<std::endl;
	__n = mat.cols();
	SAFE_DELETE_ARRAY(__d);
	SAFE_DELETE_ARRAY(__e);
	SAFE_DELETE_ARRAY(__tau);
	__d = new scalar[__n];
	__e = new scalar[__n-1];
	__tau = new scalar[__n-1];
	copt_lapack_sytrd(__uplo,__n,const_cast<scalar*>(mat.dataPtr()),mat.rows(),__d,__e,__tau,&__info);
}

template<class Matrix>
SymmetricToTridiagonal<Matrix>::~SymmetricToTridiagonal()
{
	SAFE_DELETE_ARRAY(__d);
	SAFE_DELETE_ARRAY(__e);
}

template<class Matrix>
typename SymmetricToTridiagonal<Matrix>::index SymmetricToTridiagonal<Matrix>::n() const
{
	return __n;
}

template<class Matrix>
const typename SymmetricToTridiagonal<Matrix>::scalar* SymmetricToTridiagonal<Matrix>::d() const
{
	return __d;
}

template<class Matrix>
const typename SymmetricToTridiagonal<Matrix>::scalar* SymmetricToTridiagonal<Matrix>::e() const
{
	return __e;
}

/*********************Implementation of 'PartialEigenSolver'****************/

template<class Matrix>
PartialEigenSolver<Matrix>::PartialEigenSolver()
	:
	__w(NULL),
	__iblock(NULL),
	__isplit(NULL),
	__is_solved(false)
{
}

template<class Matrix>
PartialEigenSolver<Matrix>::PartialEigenSolver(const Matrix& mat)
	:
	__w(NULL),
	__iblock(NULL),
	__isplit(NULL),
	__is_solved(false)
{
	compute(mat);
}

template<class Matrix>
PartialEigenSolver<Matrix>::~PartialEigenSolver()
{
	SAFE_DELETE_ARRAY(__w);
	SAFE_DELETE_ARRAY(__iblock);
	SAFE_DELETE_ARRAY(__isplit);
}

template<class Matrix>
void PartialEigenSolver<Matrix>::compute( const Matrix& mat )
{
	if (!mat.isSymmetric())
		std::cerr<<" SymmetricToTridiagonal Warning: Please make sure that whether the input matrix is symmetric or not!"<<std::endl;
	__tri.compute(mat);
	__dim = mat.cols();
	SAFE_DELETE_ARRAY(__w);
	SAFE_DELETE_ARRAY(__iblock);
	SAFE_DELETE_ARRAY(__isplit);
	__w = new scalar[__dim];
	__iblock = new index[__dim];
	__isplit = new index[__dim];
	__is_solved = false;
}

template<class Matrix>
void PartialEigenSolver<Matrix>::solveIndex( index il , index iu)
{
	scalar *work = new scalar[4*__dim];
	index *iwork = new index[3*__dim];
	copt_lapack_stebz('I','B',__dim,0.0,0.0,il,iu,0.0,const_cast<scalar*>(__tri.d()),const_cast<scalar*>(__tri.e()) ,&__m,&__nsplit,__w,__iblock,__isplit,work,iwork,&__info);
	delete[]work;
	delete[]iwork;
	__is_solved = true;
}

template<class Matrix>
typename PartialEigenSolver<Matrix>::scalar PartialEigenSolver<Matrix>::computeLargestEigenvalue()
{
	solveIndex(__dim,__dim);
	return __w[0];
}


}
#endif