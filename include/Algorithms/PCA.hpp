#ifndef PCA_HPP__
#define PCA_HPP__

namespace COPT
{

template<class InputType,class OutputType=InputType>
class RedesignedSolver
{
private:
	virtual void doCompute(const InputType&) = 0;
	virtual OutputType doSolve(const InputType&) = 0;
public:
	virtual void compute(const InputType& input){
		this->doCompute(input);
	}

	virtual OutputType solve(const InputType& input){
		return this->doSolve(input);
	}
};

/** Principle Component Analysis(PCA), a popular tool for high-dimensional data processing
  * assumes that given high-dimensional data lies near a much lower-dimensional linear subspace.
  * The solver solves the problem:
  * 		min_{A,E} \|E\|_F subject to rank(A)<=r, D=A+E
  * However, the approach for solving the problem is quite simple. z
  */
template<class Matrix>
class PCA
	:
	public RedesignedSolver<Matrix,typename Matrix::DMatrix>
{
public:
	enum Direction{
		Row,
		Col
	};
private:

	typedef typename Matrix::scalar 			scalar;
	typedef typename Matrix::index 				index;
	typedef typename Matrix::DMatrix 			DMatrix;
	typedef typename Matrix::DVector 			DVector;

	/** the direction of the samples */
	Direction 			__dir;
	DMatrix 			__m;
	DVector 			__val;
	DMatrix 			__v;
	index 				__r;

	void doCompute(const Matrix& mat);
	DMatrix doSolve(const Matrix& mat);
public:
	PCA(const Matrix& m,const index r=1);
};

template<class Matrix>
PCA<Matrix>::PCA(const Matrix& m,const index r)
	:
	__dir(Row),
	__r(r)
{
	this->compute(m);
}

template<class Matrix>
void PCA<Matrix>::doCompute(const Matrix& mat)
{
	__m = mat;
	if(__dir == Row)
	{
		DVector __mean = mean(__m.rowBegin(),__m.rowEnd());
		std::cout<<"mean vector is "<<__mean<<std::endl;
		std::for_each(__m.rowBegin(),__m.rowEnd(),[&__mean](DVector& v){v=v-__mean;});
		std::cout<<"done"<<std::endl;
	}
	else if(__dir == Col)
	{
		DVector __mean = mean(__m.colBegin(),__m.colEnd());
		std::for_each(__m.colBegin(),__m.colEnd(),[&__mean](DVector& v){v=v-__mean;});
	}
}

template<class Matrix>
typename Matrix::DMatrix PCA<Matrix>::doSolve(const Matrix& mat)
{
	DMatrix mtm;
	__m.mtm(mtm);
	std::cout<<mtm<<std::endl;
	EigenSolver<Matrix> es(mtm);
	std::cout<<"eigen value "<<es.eigenValue();
	index dim = es.eigenValue().dimension();
	__val.resize(dim);
	for ( index i = 0; i < __r ; ++ i )
	{
		__val(dim-i-1) = es.eigenValue()(dim-i-1);
	}
	__v = es.eigenVector();
	return __v*DMatrix::diag(dim,dim,__val)*__v.transpose();
}

} // End of namespace COPT

#endif