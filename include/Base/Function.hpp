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

#ifndef FUNCTION_HPP__
#define FUNCTION_HPP__
/*
 Second version of 'Functions' classes
 */
namespace COPT
{

/*
 *	Scalar function
 *		input: scalar x
 *		output: scalar y
 *	template
 *		T stands for the float type that is used, for example double or float
 *
 */
template<class T>
class ScalarFunction
{
protected:

public:
	typedef T ScalarType;
	typedef T FT;
	ScalarFunction()
	{
	}
	// the deconstructor
	virtual ~ScalarFunction()
	{
	}
	// basic operation computing the value of the function
	virtual FT operator()(FT x) const = 0;
	// basic operation computing the differential of the function
	virtual FT diff(FT x) const
	{
		return ScalarDifferential<ScalarFunction>(*this).diff(x);
	}
	virtual FT sDiff(FT x) const
	{
		return ScalarDifferential<ScalarFunction>(*this).sDiff(x);
	}
};

/*
 *	'VectorFunction' returns a scalar value of a function defined on a vector
 *	call statement:
 *		double value = VectorFunction(vec):
 *			\param 'vec': input type Vector
 *			Returns the result of the function
 *		Vector diff = VectorFunction.gradient(vec):
 *			\param 'vec': input type Vector
 *			Returns the gradient of the function
 *		Matrix hessian = VectorFunction.hessian(vec):
 *			\param 'vec': input type Vector
 *			Returns the Hessian matrix of the function
 */
template<class VT>
class VectorFunction
{
protected:
	// the dimension of the problem
	int __dim;
public:
	typedef VT Vector;
	typedef typename Vector::ScalarType ScalarType;
	typedef MatrixBase<ScalarType> Matrix;
	typedef typename Vector::ScalarType FT;

	VectorFunction() :
			__dim(0)
	{
	}

	virtual ~VectorFunction()
	{
	}

	int dimension()
	{
		return __dim;
	}

	virtual FT operator()(const Vector& vec) const = 0;

	virtual Vector gradient(const Vector& vec) const
	{
		return VectorDifferential<VectorFunction>(*this).gradient(vec);
	}

	virtual Matrix hessian(const Vector& vec) const
	{
		return VectorDifferential<VectorFunction>(*this).hessian(vec);
	}
};

template<class InputArg, class OutputArg>
class RFunction
{
public:
	typedef std::function<OutputArg(const InputArg&)> FuncType;
private:
	FuncType __func;
public:

	RFunction() :
			__func([](const InputArg&)
			{	return OutputArg();})
	{
	}

	RFunction(FuncType f) :
			__func(f)
	{
	}

	virtual ~RFunction()
	{
	}

	virtual OutputArg operator ()(const InputArg& arg)
	{
		if (__func)
		{
			return __func(arg);
		}
		else
		{
			return OutputArg();
		}
	}

	RFunction& operator=(FuncType f)
	{
		__func = f;
		return *this;
	}

	RFunction& operator=(const RFunction& func)
	{
		__func = FuncType(func);
		return *this;
	}

	// apply function, the same as map function
	template<class InputIterator, class OutputIterator>
	void apply(InputIterator ib, InputIterator ie, OutputIterator ob)
	{
		std::transform(ib, ie, ob, __func);
	}

	template<class InputContainer, class OutputContainer>
	void apply(const InputContainer& input, OutputContainer& output)
	{
		output.resize(input.size());
		this->apply(input.begin(), input.end(), output.begin());
	}

	template<class InputContainer>
	VectorBase<OutputArg> apply(const InputContainer& input)
	{
		VectorBase<OutputArg> result(input.size());
		this->apply(input.begin(),input.end(),result.begin());
		return result;
	}

//	template<class InputContainer,class OutputContainer>
	auto asApplyFunction()
	{
		return ([this](auto x){return this->apply(x);});
	}

	// no class-in filter function, one can refer to a third-party filter function
};

template<class Iterator>
void filterFunc(const RFunction<typename Iterator::value_type, bool>& f,
		Iterator b, Iterator e)
{
	std::remove_if(b, e, [&f](auto x)
	{	return !f(x);});
}

template<class InputIterator, class OutputIterator>
void applyFunc(
		RFunction<typename InputIterator::value_type,
				typename OutputIterator::value_type> f, InputIterator ib,
		InputIterator ie, OutputIterator ob)
{
	std::transform(ib, ie, ob, f);
}

template<class Iterator>
typename Iterator::vaule_type reduceFunc(
		const RFunction<typename Iterator::value_type,
				typename Iterator::value_type>& f, Iterator b, Iterator e)
{
	// no initial data is given
	typename Iterator::value_type i;
	std::accumulate(b, e, i, f);
	return i;
}

} // End of namespace COPT

#endif
