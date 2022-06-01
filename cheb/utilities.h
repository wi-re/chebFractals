#pragma once
#ifndef _INCLUDE_CHEB_HEADERS_INTERNAL
#error This file is an internal header and should not be included directly. Please include <cheb/cheb.h> instead.
#endif
#include <iostream>
#include <cheb/constants.h>
#include <cheb/Interval.h>
#include <cheb/Domain.h>


namespace cheb {
	//forward declarations
	class IntervalFunction;

	namespace internal {
		// This function is used to remove any elements that are within the given
		// tolerance relative to the next element in the given vector. For example
		// 1. 1.1 2 1.1 0 with a tolerance of .2 will result in 1 2 1.1 0
		// The vector tolerance variant expects n - 1 tolerance entries.
		// The last element is always contained in the output (it has next elem).
		svec merge_duplicates(const svec& arr, const svec& tols);
		svec merge_duplicates(const svec& arr, scalar tol);
		// Given a list of intervals, this function returns a list that corresponds
		// to the sorted intervals, i.e. [1 ,0] would indicate that intervals was
		// reverse sorted beforehand. This function will also check if the sorted
		// vector would describe a proper domain, i.e. there is no gap between the
		// sorted intervals and they do not overlay (these tests are exact and not
		// defined by some tolerance threshold)
		std::vector<int32_t> sortindex(const std::vector<Interval>& intervals);
	}
	// Given a set of IntervalFunctions, this method will first sort the intervals they describe.
	// This also checks if the set of functions provided describes a correct domain. Once
	// the intervals are sorted the functions will be sorted accordingly and the sorted
	// list of functions will be returned.
	std::vector<IntervalFunction> check_funs(const std::vector<IntervalFunction>& funs);
	// Given a set of IntervalFunctions, this function will evaluate all functions on their
	// respective start and end points, i.e. their support. These values are returned as
	// a <key, value> pair of arrays, i.e. the keys describe the respective x values at
	// which the funtions were evaluated to yield the returned values. This function 
	// requires the functions given to be correctly sorted, i.e. usingcheck_funs.
	std::pair<svec, svec> compute_breakdata(const std::vector<IntervalFunction>& funs);

	// Forward declarations to allow for operators at global scope for easier usage
	class IntervalFunction;
	class Function;

	IntervalFunction pow(const IntervalFunction& lhs, const IntervalFunction& rhs);
	IntervalFunction pow(const IntervalFunction& lhs, const scalar& rhs);
	IntervalFunction pow(const scalar& lhs, const IntervalFunction& rhs);

	Function pow(const Function& lhs, const Function& rhs);
	Function pow(const Function& lhs, const scalar& rhs);
	Function pow(const scalar& lhs, const Function& rhs);
}

// binary operators for IntervalFunctions and Functions, forward declared here to become friends later
cheb::IntervalFunction operator+(const cheb::IntervalFunction& lhs, const cheb::IntervalFunction& rhs);
cheb::IntervalFunction operator+(const cheb::IntervalFunction& lhs, const cheb::scalar& rhs);
cheb::IntervalFunction operator+(const cheb::scalar& rhs, const cheb::IntervalFunction& lhs);
cheb::IntervalFunction operator-(const cheb::IntervalFunction& lhs, const cheb::IntervalFunction& rhs);
cheb::IntervalFunction operator-(const cheb::IntervalFunction& lhs, const cheb::scalar& rhs);
cheb::IntervalFunction operator-(const cheb::scalar& lhs, const cheb::IntervalFunction& rhs);

cheb::IntervalFunction operator*(const cheb::IntervalFunction& lhs, const cheb::IntervalFunction& rhs);
cheb::IntervalFunction operator*(const cheb::IntervalFunction& lhs, const cheb::scalar& rhs);
cheb::IntervalFunction operator*(const cheb::scalar& lhs, const cheb::IntervalFunction& rhs);
cheb::IntervalFunction operator/(const cheb::IntervalFunction& lhs, const cheb::IntervalFunction& rhs);
cheb::IntervalFunction operator/(const cheb::IntervalFunction& lhs, const cheb::scalar& rhs);
cheb::IntervalFunction operator/(const cheb::scalar& lhs, const cheb::IntervalFunction& rhs);

cheb::Function operator+(const cheb::Function& lhs, const cheb::Function& rhs);
cheb::Function operator+(const cheb::Function& lhs, const cheb::scalar& rhs);
cheb::Function operator+(const cheb::scalar& rhs, const cheb::Function& lhs);
cheb::Function operator-(const cheb::Function& lhs, const cheb::Function& rhs);
cheb::Function operator-(const cheb::Function& lhs, const cheb::scalar& rhs);
cheb::Function operator-(const cheb::scalar& lhs, const cheb::Function& rhs);

cheb::Function operator*(const cheb::Function& lhs, const cheb::Function& rhs);
cheb::Function operator*(const cheb::Function& lhs, const cheb::scalar& rhs);
cheb::Function operator*(const cheb::scalar& lhs, const cheb::Function& rhs);
cheb::Function operator/(const cheb::Function& lhs, const cheb::Function& rhs);
cheb::Function operator/(const cheb::Function& lhs, const cheb::scalar& rhs);
cheb::Function operator/(const cheb::scalar& lhs, const cheb::Function& rhs);
