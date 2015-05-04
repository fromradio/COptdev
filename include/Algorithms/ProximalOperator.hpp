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

#ifndef PROXIMAL_OPERATORS_HPP__
#define PROXIMAL_OPERATORS_HPP__

/** 	Fast computation of famous proximal operators 
 */

namespace COPT {

/*		calculation of proximal operators
 */
template<class Function>
class ProximalOperator {

};

class QuardraticProximal {
private:
public:

};

class LogFunction {
public:
	typedef log_scalar_function_tag function_category;
	template<class FT>
	FT operator()(const FT x) {
		return -log(x);
	}
};

class AbsFunction {
public:
	typedef abs_scalar_function_tag function_category;
	template<class FT>
	FT operator()(const FT x) {
		return std::abs(x);
	}
};

template<class scalar>
scalar computeProximal(const LogFunction& func, const scalar v,
		const scalar lambda, const log_scalar_function_tag&) {
	return (v + std::sqrt(v * v + 4 * lambda)) / 2;
}

template<class scalar>
scalar computeProximal(const AbsFunction& func, const scalar v,
		const scalar lambda, const abs_scalar_function_tag&) {
	if (v >= lambda)
		return (v - lambda);
	else if (v <= -lambda)
		return (v + lambda);
	else
		return 0;
}

template<class scalar>
std::complex<scalar> computeProximal(const AbsFunction& func,
		const std::complex<scalar> v, const scalar lambda,
		const abs_scalar_function_tag&) {
	// proximal operator on real and imag part
	scalar r, i;
	if (v.real() >= lambda)
		r = v.real() - lambda;
	else if (v.real() <= -lambda)
		r = v.real() + lambda;
	else
		r = 0;
	if (v.imag() >= lambda)
		i = v.imag() - lambda;
	else if (v.imag() <= -lambda)
		i = v.imag() + lambda;
	else
		i = 0;
	return std::complex<scalar>(r, i);
}

template<class Vector, class scalar>
void computeProximal(const AbsFunction& func, const Vector& v,
		const scalar lambda, const abs_scalar_function_tag&, Vector& x) {
	if (x.size() != v.size())
		x.resize(v.size());
	for (int i = 0; i < v.size(); ++i) {
		x(i) = computeProximal(func, v(i), lambda);
	}
}

template<class Function, class scalar>
scalar computeProximal(const Function& func, const scalar v,
		const scalar lambda) {
	return computeProximal(func, v, lambda,
			typename Function::function_category());
}

template<class Function, class scalar>
std::complex<scalar> computeProximal(const Function& func,
		const std::complex<scalar>& v, const scalar lambda) {
	return computeProximal(func, v, lambda,
			typename Function::function_category());
}

template<class Function, class Vector, class scalar>
void computeProximal(const Function& func, const Vector& v, const scalar lambda,
		Vector& x) {
	computeProximal(func, v, lambda, typename Function::function_category(), x);
}
}

#endif
