#include <cheb/cheb.h>
#include <cheb/utilities.h>
#include <cheb/Function.h>

namespace cheb {
	std::size_t Domain::size() const {
		return breakpoints.size();
	}
	Domain::Domain(svec _breakpoints) :breakpoints(_breakpoints) {
		if (breakpoints.size() == 0) return;
		bool invalidRegion = breakpoints.size() == 1;
		if (breakpoints.size() > 0)
			for (int32_t i = 0; i < breakpoints.size() - 1; ++i)
				invalidRegion = invalidRegion || breakpoints[i + 1] - breakpoints[i] <= 0;
		if (invalidRegion)
			throw std::invalid_argument("Domain with gaps provided");
	}
	std::vector<Interval> Domain::intervals() const {
		std::vector<Interval> intervals;
		for (int32_t i = 0; i < breakpoints.size() - 1; ++i)
			intervals.emplace_back(breakpoints[i], breakpoints[i + 1]);
		return intervals;
	}

	Interval Domain::support() const {
		return Interval(*std::begin(breakpoints), *(std::end(breakpoints) - 1));
	}
	bool Domain::contains(const Domain& other) const {
		auto [a, b] = support();
		auto [x, y] = other.support();
		auto bounds = vec2(1. - HTOL, 1. + HTOL);
		auto lbnd = std::min(a * bounds.first, a * bounds.second);
		auto rbnd = std::max(b * bounds.first, b * bounds.second);
		return (lbnd <= x) && (y <= rbnd);
	}

	Domain Domain::united(const Domain& other) const {
		auto [a, b] = support();
		auto [x, y] = other.support();
		vec2 dspt = vec2(abs(a - x), abs(b - y));
		vec2 htol = vec2(std::max(HTOL, HTOL * abs(a)), std::max(HTOL, HTOL * abs(b)));
		if (dspt.first > htol.first || dspt.second > htol.second) {
			throw std::invalid_argument("Supports of provided domains do no match");
		}
		return merged(other);
	}
	Domain Domain::merged(const cheb::svec& other) const {
		svec all_breakpoints = concatenate(breakpoints, other);
		svec unique = uniqueElements(all_breakpoints);
		svec mergetol = unique;
		for (auto& m : mergetol)
			m = std::max(HTOL, HTOL * std::abs(m));
		//auto mergetol = evaluate(unique, [](scalar val) {return std::max(HTOL, HTOL * abs(val)); });
		auto mgd_breakpoints = internal::merge_duplicates(unique, mergetol);
		return Domain(mgd_breakpoints);
	}
	Domain Domain::merged(const Domain& other) const {
		svec all_breakpoints = concatenate(breakpoints, other.breakpoints);
		svec unique = uniqueElements(all_breakpoints);
		svec mergetol = unique;
		for (auto& m : mergetol)
			m = std::max(HTOL, HTOL * std::abs(m));
		//auto mergetol = evaluate(unique, [](scalar val) {return std::max(HTOL, HTOL * abs(val)); });
		auto mgd_breakpoints = internal::merge_duplicates(unique, mergetol);
		return Domain(mgd_breakpoints);
	}
	Domain Domain::restricted(const Domain& other) const {
		if (!contains(other))
			throw std::invalid_argument("Not Subdomain");
		auto dom = merged(other);
		auto [a, b] = other.support();
		auto bounds = vec2(1. - HTOL, 1. + HTOL);
		auto lbnd = std::min(a * bounds.first, a * bounds.second);
		auto rbnd = std::max(b * bounds.first, b * bounds.second);
		svec newDoms;
		for (auto b : dom.breakpoints)
			if (lbnd <= b && b <= rbnd)
				newDoms.push_back(b);
		return Domain(newDoms);
	}
	bool Domain::operator==(const Domain& other) const {
		if (size() != other.size())
			return false;
		auto dpdt = cheb::absolute(vSub(breakpoints,other.breakpoints));
		auto htol = std::max(HTOL, HTOL * std::abs(cheb::maximum(breakpoints)));
		bool flag = true;
		for (auto d : dpdt)
			flag = flag && d <= htol;
		return flag;
	}
	bool Domain::operator!=(const Domain& other) const {
		return !(*this == other);
	}
	scalar Domain::operator[](std::size_t idx) const {
		return breakpoints[idx];
	}
}