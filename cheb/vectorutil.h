#pragma once 
#ifndef _INCLUDE_CHEB_HEADERS_INTERNAL
#error This file is an internal header and should not be included directly. Please include <cheb/cheb.h> instead.
#endif
#include <vector>
#include <complex>
#include <functional>
#include <algorithm>
#include <numbers>
#include <random>
#include <sstream>
#include <iostream>
#include <cheb/constants.h>
#include <iomanip>

namespace cheb {
	// vector utility functions, i.e. conversion

	// convert one iterable container to another type
	// Uy needs to be manually specified when calling this function
	template<typename Uy, typename Ty>
	inline auto convert(const Ty& v) {
		Uy result(v.size());
		for (int32_t i = 0; i < v.size(); ++i)
			result[i] = v[i];
		return result;
	}
	// convert one iterable container to another after applying a
	// scalar function to every element
	template<typename Uy, typename Ty, typename Func>
	inline auto transform(const Ty& v, Func&& fn) {
		Uy result(v.size());
		for (int32_t i = 0; i < v.size(); ++i)
			result[i] = fn(v[i]);
		return result;
	}
	// function used to evaluate a function taking a Ty and returning a Ty on a vector of Tys
	template<typename Ty, typename Func>
	inline auto evaluate(const std::vector<Ty>& vals, Func&& fn) {
		return transform<std::vector<std::decay_t<decltype(fn(vals[0]))>>>(vals, fn);
	}
	// concatenate two iterable instances of an iterable container together
	template<typename Ty = svec>
	inline auto concatenate(const Ty& lhs, const Ty& rhs) {
		auto n = lhs.size();
		auto m = rhs.size();
		Ty result(n + m);
		for (int32_t i = 0; i < n; ++i)
			result[i] = lhs[i];
		for (int32_t i = 0; i < m; ++i)
			result[n + i] = rhs[i];
		return result;
	}
	// remove all duplicate elements from an iterable container
	// works by first convertign the input to a vector
	template<typename Ty = scalar>
	inline auto uniqueElements(const Ty& va) {
		std::vector vec(std::begin(va), std::end(va));
		std::sort(std::begin(vec), std::end(vec));
		vec.erase(std::unique(std::begin(vec), std::end(vec)), std::end(vec));
		Ty result(vec.size());
		for (int32_t i = 0; i < vec.size(); ++i)
			result[i] = vec[i];
		return result;
	}
	// zip two iterable containers together and return a vector containing pairs of both 
	// input arguments. requires size to be equal for both to work.
	template<typename Ty, typename Uy>
	inline auto zip(const Ty& lhs, const Uy& rhs) {
		assert(lhs.size() == rhs.size());
		using T = std::decay_t<decltype(lhs[0])>;
		using U = std::decay_t<decltype(rhs[0])>;
		std::vector<std::pair<T, U>> res(lhs.size());
		for (int32_t i = 0; i < res.size(); ++i)
			res[i] = std::make_pair(lhs[i], rhs[i]);
		return res;
	}
	// return the first n elements of an iterable container
	// if n < 0 the container is cut to the n-th last element 
	// i.e. n = -1 removes the last element (same as n = ak.size() -1)
	template<typename T>
	inline auto subVec(const T& ak, int32_t lastIdx) {
		using Ty = std::decay_t<decltype(ak[lastIdx])>;
		if (lastIdx >= 0)
			return convert<T>(std::vector<Ty>{ std::begin(ak), std::begin(ak) + lastIdx });
		return convert<T>(std::vector<Ty>{ std::begin(ak), std::end(ak) - 1 - lastIdx });
	}
	// return a sorted copy of the input iterable container
	template<typename Ty >
	inline auto sort(Ty va) {
		std::sort(std::begin(va), std::end(va));
		return va;
	}
	// return a reversed copy of the input container
	template<typename Ty>
	inline auto flip(const Ty& v) {
		Ty result(v.size());
		for (int32_t i = 0; i < v.size(); ++i)
			result[v.size() - 1 - i] = v[i];
		return result;
	}
	// returns the index of the smallest element in v
	template<typename Ty>
	inline auto argmin(const Ty& v) {
		using T = std::decay_t<decltype(v[0])>;
		T min = v[0];
		int32_t idx = 0;
		for (int32_t i = 0; i < v.size(); ++i) {
			if (v[i] < min) {
				idx = i;
				min = v[i];
			}
		}
		return idx;
	}
	// returns a slice of an array from the entry min to the entry max with the given stride
	// this function returns a slice of the input container (reference)
	//template<typename T = svec>
	//inline auto generateSlice(T& va, int32_t min, int32_t max = INT_MAX, int32_t stride = 1) {
	//	if (max == INT_MAX)
	//		max = (len)va.size();
	//	auto elems = (max - min) / stride;
	//	//elems = elems + (max - min % abs(stride == 0 ? 0 : 1);
	//	auto s = std::slice(min, stride < 0 ? std::max(elems, 0) : std::max(elems, 0), stride);
	//	return va[s];
	//}
	//// returns a slice of an array from the entry min to the entry max with the given stride
	//// this function returns a copy of the input container
	//template<typename T = svec>
	//inline auto slice(T& va, int32_t min, int32_t max = INT_MAX, int32_t stride = 1) {
	//	return T(generateSlice(va, min, max, stride));
	//}

	// vector math functions, i.e. scan

	// function that either check if all elements of a vbool are true
	// or apply a func transforming a scalar or complex to bool and then 
	// check if all elements are true. see also any
	bool all(const bvec& cfs);
	bool all(const svec& cfs, std::function<bool(scalar)> predicate);
	bool all(const cvec& cfs, std::function<bool(complex)> predicate);
	// equivalents to all that check if any element of a vbool is true or
	// check if any element of a vector is set to true after applying a 
	// predicate to all elements
	bool any(const bvec& cfs);
	bool any(const svec& cfs, std::function<bool(scalar)> func);
	bool any(const cvec& cfs, std::function<bool(complex)> func);
	// evaluates the L_infinite norm on a svec, i.e. the maximum absolute element
	scalar infnorm(const svec& vals);
	scalar infnorm(const svec& a, const svec& b);
	// return the euclidean product of two iterable containers (dot product)
	template<typename Ty, typename Uy>
	inline auto dot(const Ty& va, const Uy& vb) {
		using T = std::decay_t<decltype(va[0] * vb[0])>;
		T sum = 0.0;
		for (int32_t i = 0; i < va.size(); ++i)
			sum += va[i] * vb[i];
		return sum;
	}
	// perform a scan operation binop in the rovided container
	// the first entry contains binop(init, first element)!
	template<typename Ty, typename Func>
	inline auto scan(const Ty& v, Func&& binop, std::decay_t<decltype(v[0])> init = 0.0) {
		Ty result(v.size());
		result[0] = binop(init, v[0]);
		for (int32_t i = 1; i < v.size(); ++i)
			result[i] = binop(v[i], result[i - 1]);
		return result;
	}
	// evaluates the difference between consecutive elements of a vector
	// the returned vector is shorter by 1 element w.r.t the input
	template<typename Ty>
	auto diff(const Ty& vec) {
		Ty result(vec.size() - 1);
		for (int32_t i = 0; i < vec.size() - 1; ++i)
			result[i] = vec[i + 1] - vec[i];
		return result;
	}


	// vector generation functions

	// return a svec filled with random values in the range [0,1]
	// uses a static internal state based on random_device 
	svec randn(int32_t n);
	// returns a svec with elements from the inclusive range [min, max]
	// requires n to be at least 2.
	svec linspace(scalar min, scalar max, int32_t n);
	svec linspaceV(scalar min, scalar max, int32_t n);
	// returns a lvec containing values from the 
	// inclusive range [min, max] in increments of 1
	lvec arange(len min, len max);

	// utility functions

	// converts a given scalar function to a vector function
	vfunc funcToVfunc(func fn);
	// prints an iterable container on cout as [ elems... ]
	template<typename Ty>
	auto printVec(const Ty& vec) {
		std::stringstream sstream;
		std::cout << "[ ";
		if (vec.size() > 0)
			for (int32_t i = 0; i < vec.size() - 1; ++i)
				std::cout << std::setprecision(17) << vec[i] << " ";
		if (vec.size() > 0)
			std::cout << std::setprecision(17) << (*(std::end(vec) - 1));
		std::cout << "]\n";
	}

#define OPERATION(name, op)                                             \
  template <typename T, typename U>                                     \
  auto v##name(const std::vector<T>& lhs, const std::vector<U>& rhs) {  \
    if (lhs.size() != rhs.size())                                       \
      throw std::invalid_argument("Length mismatch");                   \
    using V = decltype(lhs[0] op rhs[0]);                               \
    std::vector<V> res(lhs.size(), V());                                \
    for (int32_t i = 0; i < lhs.size(); ++i) res[i] = lhs[i] op rhs[i]; \
    return res;                                                         \
  }                                                                     \
  template <typename T>                                     \
  auto v##name(const T& lhs, const std::vector<T>& rhs) {               \
    using V = decltype(lhs op rhs[0]);                                  \
    std::vector<V> res(rhs.size(), V());                                \
    for (int32_t i = 0; i < rhs.size(); ++i) res[i] = lhs op rhs[i];    \
    return res;                                                         \
  }                                                                     \
  template <typename T>                                     \
  auto v##name(const std::vector<T>& lhs, const T& rhs) {               \
    using V = decltype(lhs[0] op rhs);                                  \
    std::vector<V> res(lhs.size(), V());                                \
    for (int32_t i = 0; i < lhs.size(); ++i) res[i] = lhs[i] op rhs;    \
    return res;                                                         \
  }

OPERATION(Multiply, *);
OPERATION(Divide, /);
OPERATION(Add, +);
OPERATION(Sub, -);

template<typename T>
auto cos(const std::vector<T>& arg) {
	auto r = arg;
	for (auto& v : r)
		v = std::cos(v);
	return r;
}
template<typename T>
auto sin(const std::vector<T>& arg) {
	auto r = arg;
	for (auto& v : r)
		v = std::sin(v);
	return r;
}
template<typename T>
auto log10(const std::vector<T>& arg) {
	auto r = arg;
	for (auto& v : r)
		v = std::log10(v);
	return r;
}

template<typename T>
auto reverse(std::vector<T>& arg) {
	for (int32_t i = 0; i < arg.size() / 2; ++i)
		std::swap(arg[i], arg[arg.size() - 1 - i]);
}

template<typename T>
auto sum(std::vector<T> arg) {
	T sum = T(0);
	for (auto v : arg)
		sum += v;
	return sum;
}
template<typename T>
auto absolute(const std::vector<T>& arg) {
	auto r = arg;
	for (auto& v : r)
		v = std::abs(v);
	return r;
}

template<typename T>
auto maximum(const std::vector<T>& arg) {
	T max = -std::numeric_limits<T>::max();
	for (auto v : arg)
		max = std::max(v,max);
	return max;
}
template<typename T>
auto power(const std::vector<T>& arg, T power) {
	auto r = arg;
	for (auto& v : r)
		v = std::pow(v, power);
	return r;
}
template<typename T>
auto power(T arg, std::vector<T> power) {
	auto r = power;
	for (auto& v : r)
		v = std::pow(arg, v);
	return r;
}
template<typename T>
auto power(std::vector<T> arg, std::vector<T> power) {
	auto r = power;
	for(int32_t i = 0; i < power.size(); ++i)
		r[i] = std::pow(arg[i], power[i]);
	return r;
}


	// helper macros for debugging
#define debugPrintVec(x) printf("%d - %s [%s] = ", __LINE__, #x, typeid(decltype(x)).name());cheb::printVec(x);
#define debugPrint(x) printf("%d - %s [%s] = %s\n", __LINE__, #x, typeid(decltype(x)).name(),std::to_string(x).c_str());
}
