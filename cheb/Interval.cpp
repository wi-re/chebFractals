#include <cheb/cheb.h>
#include <cheb/utilities.h>
#include <cheb/Function.h>

namespace cheb {
	bool Interval::contains(const Interval& other) const {
		auto [x, y] = other;
		return (a <= x) && (y <= b);
	}
	bvec Interval::isinterior(const svec& x) const {
		bvec res(x.size());
		for(int32_t i = 0; i < x.size(); ++i)
			res[i] = (a < x[i]) && (x[i] < b);
		return res;
	}
	bool operator==(const Interval& lhs, const Interval& rhs) {
		auto [a, b] = lhs;
		auto [x, y] = rhs;
		return (a == x) && (b == y);
	}
	bool operator!=(const Interval& lhs, const Interval& rhs) {
		return !(lhs == rhs);
	}
}