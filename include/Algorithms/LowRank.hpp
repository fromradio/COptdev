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


#ifndef LOW_RANK_HPP__
#define LOW_RANK_HPP__

namespace COPT
{

template<class kernel, class Time = NoTimeStatistics>
class APGSolver
    :
    public GeneralSolver<kernel, Time, typename kernel::Matrix>
{
private:
	typedef typename kernel::index		    index;
	typedef typename kernel::scalar 		scalar;
	typedef typename kernel::Vector 		Vector;
	typedef typename kernel::Matrix 		Matrix;

	const Matrix&        __D;
	const scalar         __lam;
	index                __maxiteration;

	index                __dr;
	index                __dc;

	scalar               __mu;
	scalar               __mu_bar;
	scalar               __delta;
	scalar               __eta;
	scalar               __t_forward;
	scalar               __t;
	Matrix               __A_forward;
	Matrix               __A;
	Matrix               __E_forward;
	Matrix               __E;

	APGSolver();
	void solvingBegin();
	void doCompute();
	void funcS(const scalar, Matrix&);
	scalar max(scalar, scalar);
	scalar doOneIteration();

public:
	APGSolver(
		const Matrix&    D,
		const scalar     lam,
		const index maxiteration = 10000);

	const Matrix& result() const;
	const Matrix& result_E() const;

	scalar objective() const;

};


template<class kernel, class Time>
void APGSolver<kernel, Time>::doCompute()
{
	this->__max_iteration = __maxiteration;
    __mu = 0.99*__D.operationNorm();

    // std::cout<<"original mu:"<<std::endl;
    // std::cout<<__mu<<std::endl;

    __delta = 1e-9;
    __mu_bar = __delta*__mu;
    __dr = __D.rows();
    __dc = __D.cols();
}

template<class kernel, class Time>
APGSolver<kernel, Time>::APGSolver(
	const Matrix& D,
	const scalar lam,
	const index maxiteration)
    :
    __D(D),
    __lam(lam),
    __maxiteration(maxiteration)
{
	this->doCompute();
}

template<class kernel, class Time>
void APGSolver<kernel, Time>::solvingBegin()
{
	__A_forward.resize(__dr, __dc);
	__A.resize(__dr, __dc);
	__E_forward.resize(__dr, __dc);
	__E.resize(__dr, __dc);
	__t_forward = 1;
	__t = 1;
	__eta = 0.9;
}

template<class kernel, class Time>
void APGSolver<kernel, Time>::funcS(const scalar ep, Matrix& S)
{
	index r = S.rows(), c = S.cols();
	for(int i = 0; i < r; i++)
		for(int j = 0; j < c; j++)
			if(S(i,j) > ep)
			    S(i,j) = S(i,j) - ep;
			else if(S(i,j) < -ep)
				S(i,j) = S(i,j) + ep;
			else
				S(i,j) = 0;
}

template<class kernel, class Time>
typename APGSolver<kernel, Time>::scalar APGSolver<kernel, Time>::max(scalar x, scalar y)
{
	if(x >= y)
		return x;
	else
		return y;
}

template<class kernel, class Time>
typename APGSolver<kernel, Time>::scalar APGSolver<kernel, Time>::doOneIteration()
{
	Matrix               __YA;
	Matrix               __YE;
	Matrix               __GA;
	Matrix               __GE;
	Matrix               __SA;
	Matrix               __SE;
	scalar               s;

	Matrix               tempS; 

	__YA = __A + (__t_forward - 1)/__t*(__A - __A_forward);
	__YE = __E + (__t_forward - 1)/__t*(__E - __E_forward);

	__GA = __YA - 0.5*(__YA + __YE - __D);
	COPT::SVD<Matrix> svd(__GA);
	__A_forward = __A;
    tempS = svd.S();

    // std::cout<<"svd.s"<<std::endl;
    // std::cout<<tempS<<std::endl;

    funcS(__mu/2,tempS);

 //    std::cout<<"ep"<<std::endl;
	// std::cout<<__mu/2<<std::endl;
	// std::cout<<"shit!"<<std::endl;
	// std::cout<<tempS<<std::endl;

	__A = svd.U()*tempS*svd.VT();

	__GE = __YE - 0.5*(__YA + __YE - __D);
	__E_forward = __E;
	funcS(__lam*__mu/2,__GE);
	__E = __GE;

	__t_forward = __t;
	__t = (1 + sqrt(1 + 4*__t*__t))/2;

    // std::cout<<"why?"<<std::endl;
    // std::cout<<__eta*__mu<<std::endl;
    // std::cout<<__mu_bar<<std::endl;

	__mu = max(__eta*__mu, __mu_bar);


	// std::cout<<"fucking mu?"<<std::endl;
	// std::cout<<__mu<<std::endl;

	__SA = 2*(__YA - __A) + (__A + __E - __YA - __YE);
	__SE = 2*(__YE - __E) + (__A + __E - __YA - __YE);
	s = sqrt(pow(__SA.frobeniusNorm(),2) + pow(__SE.frobeniusNorm(),2));

	// std::cout<<"thresh?"<<std::endl;
	// std::cout<<this->__thresh<<std::endl;
	// std::cout<<"WTF!"<<std::endl;
	// std::cout<<s<<std::endl;

	return s;
}

template<class kernel, class Time>
const typename APGSolver<kernel, Time>::Matrix& APGSolver<kernel, Time>::result() const
{
	return __A;
}

template<class kernel, class Time>
const typename APGSolver<kernel, Time>::Matrix& APGSolver<kernel, Time>::result_E() const
{
	return __E;
}

template<class kernel, class Time>
typename APGSolver<kernel, Time>::scalar APGSolver<kernel, Time>::objective() const
{
	// return (__A + __E - __D).operationNorm();
	return 0;
}



template<class InputType,class OutputType = InputType,class ConjugateType = InputType>
class COPTMap
	:
	public COPTObject
{
public:
	virtual ~COPTMap(){}

	virtual OutputType operator() (const InputType& ) = 0;

	virtual ConjugateType conjugate( const OutputType& ) = 0;
};

template<class Matrix>
class CompletionMap
	:
	public COPTMap<Matrix,typename Matrix::DVector>
{
private:
	typedef typename Matrix::scalar 		scalar;
	typedef typename Matrix::index 			index;
	typedef typename Matrix::DVector 		DVector;
	typedef COPT::SpMatrixBase<scalar,index>	SpMatrix;

	index 				__nr;
	index 				__nc;

	index 				__nnz;
	index 				*__ir;
	index 				*__jc;

public:

	CompletionMap() = delete;
	CompletionMap(const index nr, const index nc, const index nnz, const index *ir, const index *jc);
	CompletionMap(const CompletionMap& cm) = delete;

	DVector operator() (const Matrix& mat);
	Matrix conjugate (const DVector& vec);

};

template<class Matrix>
CompletionMap<Matrix>::CompletionMap(const index nr, const index nc, const index nnz, const index *ir, const index *jc)
	:
	__nr(nr),
	__nc(nc),
	__nnz(nnz),
	__ir(nullptr),
	__jc(nullptr)
{
	__ir = new index[nnz];
	blas::copt_blas_copy(nnz,ir,1,__ir,1);
	__jc = new index[nc+1];
	blas::copt_blas_copy(nc+1,jc,1,__jc,1);
}

template<class Matrix>
typename CompletionMap<Matrix>::DVector CompletionMap<Matrix>::operator() (const Matrix& mat)
{
	DVector result(__nnz);
	index i;
	for(index c = 0 ; c < __nc ; ++ c )
	{
		for (i = __jc[c]; i < __jc[c+1]; ++ i ){
			result(i) = mat(__ir[i],c);
		}
	}
	return result;
}

template<class Matrix>
Matrix CompletionMap<Matrix>::conjugate(const DVector& vec)
{
	assert(vec.size()==__nnz);
	return SpMatrix(__nr,__nc,__nnz,__jc,__ir,vec.dataPtr()).toDenseMatrix();
}

template<class kernel>
class APGLSolver
	:
	public GeneralSolver<kernel,NoTimeStatistics,typename kernel::Matrix>
{
private:

	/** 
	**/
	typedef typename kernel::scalar 			scalar;
	typedef typename kernel::podscalar 			podscalar;
	typedef typename kernel::index 				index;
	typedef typename kernel::Matrix 			Matrix;
	typedef typename kernel::Vector 			Vector;
	typedef typename kernel::SpMatrix 			SpMatrix;

	index 						__m;
	index 						__n;
	/** completion map */
	CompletionMap<Matrix> 		__cm;
	Vector 						__b;
	/** beta */
	scalar 						__beta;
	/** mu */
	scalar 						__mu;
	/** x_{k-1} */
	Matrix 						__x_m;
	/** x_k */
	Matrix 						__x;
	/** t_{k-1} */
	scalar 						__t_m;
	/** t_k */
	scalar 						__t;
	/** tau */
	scalar 						__tau;
	/** objective function */
	scalar 						__o;

	void solvingBegin();
	void doCompute();
	podscalar doOneIteration();

public:

	APGLSolver(const index nr, const index nc, const index nnz, const index *ir, const index *jc, const Vector& b, const scalar beta = 0.5, const scalar mu = 0.005);
	APGLSolver(const SpMatrix& spmat, const scalar beta = 2.0, const scalar mu = 0.000005);

	podscalar objective() const;

	const Matrix& result() const;

	
};

template<class kernel>
APGLSolver<kernel>::APGLSolver(
	const index nr, const index nc, const index nnz,
	const index *ir, const index *jc, const Vector& b,
	const scalar beta, const scalar mu)
	:

	__m(nr),
	__n(nc),
	__cm(nr,nc,nnz,ir,jc),
	__b(b),
	__beta(beta),
	__mu(mu)
{
	std::cout<<"here"<<std::endl;
	assert(b.size()==nnz);
	this->compute();
}

template<class kernel>
APGLSolver<kernel>::APGLSolver(
	const SpMatrix& spmat, const scalar beta, const scalar mu)
	:
	__m(spmat.rows()),
	__n(spmat.cols()),
	__cm(spmat.rows(),spmat.cols(),spmat.elementSize(),spmat.rowIndex(),spmat.columnPointer()),
	__beta(beta),
	__mu(mu)
{
	__b.setArray(spmat.elementSize(),spmat.values());
}

template<class kernel>
void APGLSolver<kernel>::doCompute()
{
}

template<class kernel>
void APGLSolver<kernel>::solvingBegin()
{
	__x.resize(__m,__n);
	__x_m.resize(__m,__n);
	__t = 1.0;
	__t_m = 1.0;
	__tau = 2.0;
	__o = 0.0;
	std::cout<<"solving begin!"<<std::endl;
}

template<class kernel>
typename APGLSolver<kernel>::podscalar APGLSolver<kernel>::doOneIteration()
{

	Matrix y = __x+((__t_m-1)/__t)*(__x-__x_m);
	Matrix g = y - (1.0/__tau)*__cm.conjugate(__cm(y)-__b);
	SVD<Matrix> de(g);
	Matrix S = de.S();
	index mi = std::min(S.rows(),S.cols());
	for ( index i = 0 ; i < mi ; ++ i )
	{
		if(S(i,i)>__mu/__tau)
			S(i,i) = S(i,i)-__mu/__tau;
		else
			S(i,i) = 0;
	}
	__x_m = __x;
	__x = de.U()*de.S()*de.VT();
	__t_m = __t;
	__t = (1.0+std::sqrt(1+4*__t*__t))/2.0;

	scalar ma = std::max(__x.frobeniusNorm(),1.0);
	scalar e = (__x-__x_m).frobeniusNorm()/ma;
	return e;
}

template<class kernel>
typename APGLSolver<kernel>::podscalar APGLSolver<kernel>::objective() const
{
	return __o;
}

template<class kernel>
const typename APGLSolver<kernel>::Matrix& APGLSolver<kernel>::result() const
{
	return __x;
}

}

#endif