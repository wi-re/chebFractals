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

namespace cheb {
	// Type conventions:
	// scalar is the underlying floating type, set to double
	// complex is the underlying complex floating type, set to std::complex<double>
	// len is the underlying length type, to make this work nicely with some function it is a SIGNED integer type
	// cvec types are std::vector<c> [svec lvec cvec]
	// function types are used as shorthands for funcitons taking scalars or svecs and returning them [func, vfunc]
	// vec2 is a shorthand for a pair of scalars.

	using scalar = double;
	using svec = std::vector<scalar>;

	using len = int32_t;
	using lvec = std::vector<len>;

	using complex = std::complex<scalar>;
	using cvec = std::vector<complex>;

	using func = std::function<scalar(scalar)>;
	using vfunc = std::function<svec(svec)>;

	using bvec = std::vector<bool>;
	using vec2 = std::pair<scalar, scalar>;

	// based on the define UNSIGNED_LOOPS use either signed or unsigned integers for looping over arrays
#ifndef UNSIGNED_LOOPS
	using loop_t = int32_t;
#else 
	using loop_t = std::size_t;
#endif

	// pi constant for shorthand access (requires C++20)
	constexpr auto inline pi = std::numbers::pi_v<scalar>;
	// underlying machine precision for the scalar type
	auto inline eps = 1.0 * std::numeric_limits<scalar>::epsilon();
	// a slightly higher tolerance level required for certain domain functionality
	auto inline HTOL = 5. * eps;
	// magic number used when splitting intervals during root finding.
	// this value is taking verbatim from Function 
	// see the original source: https://github.com/Function/Function/blob/614af81fd9dd64765f0ad51202b0706b0742fb09/%40chebtech/roots.m
	constexpr auto inline SPLITPOINT = -0x1.3dd6ba5c9ce03p-8;// -0.004849834917525;
	// enum to select the method to evaluate chebyshev polynomials, only used when directly calling a function
	enum struct method { clenshaw, bary };
	// The default evaluation scheme used for functions
	constexpr inline method defaultEvalScheme = method::clenshaw;
}
