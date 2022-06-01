#pragma once
#ifndef _INCLUDE_CHEB_HEADERS_INTERNAL
#error This file is an internal header and should not be included directly. Please include <cheb/cheb.h> instead.
#endif
#include <iostream>
#include <cheb/constants.h>

namespace cheb {
	/*
	This class models domains as they are used for Function. These domains are defined 
	by a set of breakpoints, i.e. [-1,0,1], and describe a set of intervals, i.e.
	the breakpoints [-1,0,1] define the intervals [-1,0] and [0,1].

	Note that the sub intervals overlap on the boundary elements. However, the isinterior
	function of intervals corrects for this aspect by refering to the exclusive interval
	when being called. Accordingly, users of this class should evaluate the respective 
	underlying function at the breakpoints separately and combine these results later.

	This class is immutable.
	*/
	class Domain {
	public:
		const svec breakpoints;
		// A domain can either be constructed from a svec describing the actual
		// breakpoints, an interval instance or using an initializer list of scalar
		// values. Domains should contain at least two breakpoints, however, it is
		// also possible to have a domain with no breakpoints and as such it is 
		// recommend to first check size() before interacting with a domain.
		Domain(svec _breakpoints = { -1.,1. });
		explicit Domain(Interval interval) :Domain(svec{ interval.a, interval.b }) {}
		Domain(std::initializer_list<scalar> list) :Domain(svec{ list }) {}

		// Construct a vector of the intervals described by the breakpoints of this
		// class. Note that these intervals are non-overlapping and consecutive
		// by constructing, i.e. they can be directly used for other functions.
		std::vector<Interval> intervals() const;
		// Returns an interval describing the lowest and highest breakpoint value
		Interval support() const;
		// Check if one domain is a superset of another, similar to intervals
		bool contains(const Domain& other) const;

		// Merges two domains together, requiring them to have equal support
		// the resulting domain will contain all breakpoints of both domains
		// without any duplicates.
		Domain merged(const cheb::svec& other) const;
		// Merges two domains together, requiring them to have equal support
		// the resulting domain will contain all breakpoints of both domains
		// without any duplicates.
		Domain merged(const Domain& other) const;
		// Merges two domains together but does not check if they have equal
		// support. This might be useful in some operations with a-prior known
		// relations of domains
		Domain united(const Domain& other) const;
		// Restricts this domain by the support of another domain, i.e. the
		// restricted domain only contains breakpoints within the support of
		// the other domain (including the other domains support)
		Domain restricted(const Domain& other) const;

		// Helper functions
		bool operator==(const Domain& other) const;
		bool operator!=(const Domain& other) const;
		std::size_t size() const;
		scalar operator[](std::size_t idx) const;
	};
}