#ifndef CONJUGATE_GRADIENT_METHOD
#define CONJUGATE_GRADIENT_METHOD

/*
 *		Powered by 'MathU'
 * 		Copyright @ MathU,
 *		Written by Ruimin Wang, ruimin.wang13@gmail.com
 */


namespace COPT{
/*
 *		A function solving linear conjugate gradient problem
 *		/param 'mat', the left hand matrix
 */
template <class Matrix,class Vector>
void conjugateGradientWithoutPrecondition(const Matrix& mat,const Vector& rhs,const Vector& dest,typename Vector::Scalar& tol){

}
}



#endif