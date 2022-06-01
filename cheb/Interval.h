#pragma once
#ifndef _INCLUDE_CHEB_HEADERS_INTERNAL
#error This file is an internal header and should not be included directly. Please include <cheb/cheb.h> instead.
#endif
#include <iostream>
#include <cheb/constants.h>

namespace cheb {
	/*
	This class provides a simple closed interval [a, b] with functions to map to
	and from this interval. The call operator maps to the interval and can be called
	with scalar value and vectors.

	Note that these intervals refer to a range from [-1,1] instead of the more common
	[0,1] as they are intended to work with CPR.
	*/
	struct Interval {
		scalar a = -1, b = 1;

		Interval(scalar _a = -1, scalar _b = 1) :a(_a), b(_b) {
			if (a >= b)
				throw std::invalid_argument("Invalid interval ( begin >= end )");
		}
		// map value from [-1., 1.] to the interval described by this instance
		inline auto mapToInterval(const scalar& y) const {
			if (a == -1. && b == 1.) return y;
			//return .5 * b * (y + 1.) + .5 * a * (1. - y);
			return std::fma(y, (.5 * (b - a)), .5 * (a + b));
		}
		inline auto mapToInterval(const svec& y) const {
			if (a == -1. && b == 1.) return y;
			auto res = y;
			for(auto& y : res)
				y = std::fma(y, (.5 * (b - a)), .5 * (a + b));
			return res;
		}
		// map value from [a, b] back to [-1., 1.]
		inline auto mapFromInterval(const scalar& x) const {
			if (a == -1. && b == 1.) return x;
			//return (2. * x - a - b) / (b - a);
			return std::fma(x, (2. / (b - a)), -(a + b) / (b - a));
		}
		inline auto mapFromInterval(const svec& x) const {
			if (a == -1. && b == 1.) return x;
			auto res = x;
			for (auto& x : res)
				x = std::fma(x, (2. / (b - a)), -(a + b) / (b - a));
			return res;
		}
		// call operator that maps from [-1.,1.] to [a, b]
		template<typename T>
		inline auto operator()(const T& x) const {
			return mapToInterval(x);
		}
		// checks if the other interval is contained in this interval
		// i.e. [-1,1] contains [-0.5,0.5]
		bool contains(const Interval& other) const;
		// returns a mask corresponding to which elements in the argument are
		// within the excolusive interval, i.e. ( a < x ) && ( x < b )
		bvec isinterior(const svec& x) const;
	};
	// helper comparison operators for intervals
	bool operator==(const Interval& lhs, const Interval& rhs);
	bool operator!=(const Interval& lhs, const Interval& rhs);
}