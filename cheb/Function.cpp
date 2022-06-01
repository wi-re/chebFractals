// #pragma once
#include <cheb/cheb.h>
#include <cheb/utilities.h>
#include <cheb/IntervalFunction.h>
#include <optional>
#include <cfloat>
#include <cstring>
#include <variant>
namespace cheb {
	Function::Function(scalar c, Domain d) {
		for (auto i : d.intervals())
			funs.emplace_back(c, i);
		funs = check_funs(funs);
		breakdata = compute_breakdata(funs);
	}
	Function::Function(std::string s, Domain d) {
		for (auto i : d.intervals())
			funs.emplace_back(i);
		funs = check_funs(funs);
		breakdata = compute_breakdata(funs);
	}
	Function::Function(Domain d) {
		for (auto i : d.intervals())
			funs.emplace_back(i);
		funs = check_funs(funs);
		breakdata = compute_breakdata(funs);
	}
	Function::Function(vfunc f, std::variant<lvec, len> n, Domain d) {
		try {
			auto domain = d;
			auto nl = std::get<lvec>(n);
			auto is = domain.intervals();
			if (nl.size() == 1) {
				for (auto i : d.intervals())
					funs.emplace_back(f, nl[0], i);
				funs = check_funs(funs);
				breakdata = compute_breakdata(funs);
				return;
			}
			if (nl.size() != is.size())
				throw std::invalid_argument("Length of polynomial lengths and length of intervals don't match\n");
			for (int32_t i = 0; i < is.size(); ++i) {
				funs.push_back(IntervalFunction(f, nl[i], is[i]));
			}
			funs = check_funs(funs);
			breakdata = compute_breakdata(funs);
		}
		catch (std::bad_variant_access&) {
			auto nl = std::get<len>(n);
			for (auto i : d.intervals())
				funs.emplace_back(f, nl, i);
			funs = check_funs(funs);
			breakdata = compute_breakdata(funs);
		}
	}
	Function::Function(vfunc f, Domain d) {
		for (auto i : d.intervals()) {
			//std::cout << "\tInterval " << i.a << " : " << i.b << std::endl;
			funs.emplace_back(f, i);
		}
		funs = check_funs(funs);
		breakdata = compute_breakdata(funs);
	}
	Function::Function(func fd, std::variant<lvec, len> n, Domain d) {
		auto f = funcToVfunc(fd);
		try {
			auto domain = d;
			auto nl = std::get<lvec>(n);
			auto is = domain.intervals();
			if (nl.size() == 1) {
				for (auto i : d.intervals())
					funs.emplace_back(f, nl[0], i);
				funs = check_funs(funs);
				breakdata = compute_breakdata(funs);
				return;
			}
			if (nl.size() != is.size())
				throw std::invalid_argument("Length of polynomial lengths and length of intervals don't match\n");
			for (int32_t i = 0; i < is.size(); ++i) {
				funs.push_back(IntervalFunction(f, nl[i], is[i]));
			}
			funs = check_funs(funs);
			breakdata = compute_breakdata(funs);
		}
		catch (std::bad_variant_access&) {
			auto nl = std::get<len>(n);
			for (auto i : d.intervals())
				funs.emplace_back(f, nl, i);
			funs = check_funs(funs);
			breakdata = compute_breakdata(funs);
		}
	}
	Function::Function(func fd, Domain d) {
		auto f = funcToVfunc(fd);
		//std::clog << "Finding polynomial approximation for domain " << d.support().a << " : " << d.support().b << std::endl;
		for (auto i : d.intervals()) {
			//std::cout << "\tInterval " << i.a << " : " << i.b << std::endl;
			//std::clog << "\tTrying to find polynomial..." << std::endl;
			std::vector<std::pair<bool, Interval>> intervals;
			std::vector<std::pair<bool, Interval>> intervalsParsed;
			intervals.push_back(std::make_pair(false, i));
			while (true) {
				intervalsParsed.clear();
				for(auto [flag, interval] : intervals)
				try {
					if (flag) {
						intervalsParsed.push_back(std::make_pair(true, interval));
							continue;
					}
					auto ifunc = IntervalFunction(f, interval);
					//std::clog << "\tFound function with degree " << ifunc.coeffs().size() << " for interval " << interval.a << " : " << interval.b << std::endl;
					funs.emplace_back(ifunc);
					intervalsParsed.push_back(std::make_pair(true, interval));
				}
				catch (std::invalid_argument iarg) {
					//std::clog << "\tCould not find function for interval " << i.a << " : " << i.b << std::endl;
					auto split = (interval.b + interval.a) * .5;// -(interval.b - interval.a) * .5 * SPLITPOINT;
					Interval i0{ interval.a,  split };
					Interval i1{ split,  interval.b };
					//std::clog << "\tSubdividing interval into " << i0.a << " : " << i0.b << " x " << i1.a << " : " << i1.b << " -> "  << i0.b - i0.a << " : "  << i1.b - i1.a << std::endl;
					intervalsParsed.push_back(std::make_pair(false, i0));
					intervalsParsed.push_back(std::make_pair(false, i1));
				}
				std::swap(intervals, intervalsParsed);
				bool done = true;
				for (auto is : intervals) {
					done = done && is.first;
				}
				if (done) break;
			}
		}
		funs = check_funs(funs);
		breakdata = compute_breakdata(funs);
	}
	Function::Function(std::vector<IntervalFunction> _funs) {
		funs = check_funs(_funs);
		breakdata = compute_breakdata(funs);
	}
	Function::Function() {
		funs = check_funs({});
		breakdata = compute_breakdata(funs);
	}

	svec Function::operator()(svec x, method how) const {
		if (isempty()) return svec{};
		auto out = svec(x.size(), std::nan(""));
		// Interior point values
		for (const auto& fun : funs) {
			auto idx = fun.interval().isinterior(x);
			for (int32_t i = 0; i < x.size(); ++i)
				out[i] = idx[i] ? fun(x[i], how) : out[i];
		}
		// Brekapoint values
		for (int32_t b = 0; b < breakdata.first.size(); ++b) {
			for (int32_t i = 0; i < out.size(); ++i)
				out[i] = x[i] == breakdata.first[b] ? breakdata.second[b] : out[i];
		};
		// Endpoint values

		for (int32_t i = 0; i < out.size(); ++i)
			out[i] = x[i] < breakdata.first[0] ? funs[0](x[i]) : out[i];
		for (int32_t i = 0; i < out.size(); ++i)
			out[i] = x[i] > * (std::end(breakdata.first) - 1) ? (*std::rbegin(funs))(x[i]) : out[i];

		//auto lpts = x < breakdata.first[0];
		//auto rpts = x > * (std::end(breakdata.first) - 1);
		//out[lpts] = funs[0](x[lpts]);
		//out[rpts] = (*std::rbegin(funs))(x[rpts]);
		return out;
	}
	scalar Function::operator()(scalar x, method how) const {
		if (isempty()) return 0.;
		return this->operator()(svec{ x }, how)[0];
	}
	bool Function::isempty() const {
		return funs.size() == 0;
	}
	svec Function::breakpoints()  const {
		return breakdata.first;
	}
	Domain Function::domain() const {
		if (isempty()) return Domain{ svec{} };
		return Domain((breakpoints()));
	}
	Domain Function::support()const {
		if (isempty()) return Domain{ svec{} };
		return Domain{ {*std::begin(breakdata.first), *(std::end(breakdata.first) - 1)} };
	}
	scalar Function::hscale() const {
		if (isempty()) return 0.;
		return std::max(std::abs(*std::begin(breakdata.first)), std::abs(*(std::end(breakdata.first) - 1)));
	}
	bool Function::isconst() const {
		if (isempty()) return false;
		auto c = funs[0].coeffs()[0];
		for (const auto& fun : funs)
			if (!fun.isconst() || fun.coeffs()[0] != c) return false;
		return true;
	}
	scalar Function::vscale() const {
		if (isempty()) return 0.;
		scalar max = -DBL_MAX;
		for (const auto& fun : funs)
			max = std::max(max, fun.vscale());
		return max;
	}
	std::size_t Function::size() const {
		return funs.size();
	}
	Function Function::breakWith(Domain targetdomain) const {
		std::vector<IntervalFunction> newfuns;
		auto subintervals = targetdomain.intervals();
		auto intervalit = std::begin(subintervals);
		for (const auto& fun : funs)
			while (intervalit != std::end(subintervals) && fun.interval().contains(*intervalit))
				newfuns.push_back(fun.restricted(*intervalit++));
		return newfuns;
	}
	Function Function::simplified() const {
		if (isempty()) return Function();
		auto funs2 = funs;
		for (auto& fun : funs2)
			fun = fun.simplified();
		return Function(funs2);
	}
	Function Function::restricted(Domain subinterval, bool simplify) const {
		if (isempty()) return Function();
		auto newdom = domain().restricted(subinterval);
		Function fn = breakWith(newdom);;
		if (simplify) return fn.simplified();
		return fn;
	}
	svec Function::roots(bool polish, int32_t splitDegree) const {
		if (isempty()) return svec{};
		svec allrts;
		svec prvrts;
		auto htol = 1e2 * hscale() * eps;
		for (const auto& fun : funs) {
			auto rtsv = fun.roots(polish,splitDegree);
			auto rts = (rtsv);
			if (prvrts.size() > 0 && rts.size() > 0)
				if (std::abs(*std::rbegin(prvrts) - *std::begin(rts)) <= htol)
					rts = rts.size() == 1 ? svec{} : svec(std::begin(rts) + 1, std::end(rts));

			if (rts.size() > 0)
				allrts.insert(std::end(allrts), std::begin(rts), std::end(rts));
			prvrts = rts;
		}
		return (allrts);
	}
	Function Function::antiDerivative() const {
		std::vector<IntervalFunction> newfuns;
		std::optional<IntervalFunction> prevfun = std::nullopt;
		for (const auto& fun : funs) {
			auto integral = fun.antiDerivative();
			if (prevfun) {
				auto f = prevfun.value().endvalues();
				integral = integral + f[1];
			}
			newfuns.push_back(integral);
			prevfun = integral;
		}
		return Function(newfuns);
	}
	Function Function::derivative() const {
		std::vector<IntervalFunction> newfuns;
		for (const auto& fun : funs)
			newfuns.push_back(fun.derivative());
		return Function(newfuns);
	}
	scalar Function::definiteIntegral() const {
		scalar sum = 0.;
		for (const auto& fun : funs)
			sum += fun.defIntegral();
		return sum;
	}
	Function Function::operator-() {
		if (isempty()) return Function();
		std::vector<IntervalFunction> funs;
		for (auto& fun : this->funs)
			funs.push_back(-fun);
		return Function(funs);
	}
	scalar Function::dot(const Function& f) const {
		if (isempty()) return 0.;
		return (*this * f).definiteIntegral();
	}
	Function Function::absolute() {
		if (isempty()) return Function();
		auto olddom = domain();
		auto rts = roots();
		auto newdom = olddom.merged((rts));
		//auto newdom = domain().merged(roots());
		auto oldfuns = breakWith(newdom);
		std::vector<IntervalFunction> funs;
		for (auto& fun : oldfuns.funs)
			funs.push_back(fun.abs());
		return Function(funs);
	}
	Function Function::maximum(const Function& other) const {
		if (isempty() || other.isempty()) return Function();
		auto roots = (*this - other).roots();
		auto newdom = domain().united(other.domain()).merged((roots));
		auto interval = newdom.support();
		auto switch_ = Domain{ interval }.merged((roots));
		auto range = convert<lvec>(arange(0, (len)switch_.size() - 1));
		auto keys = vMultiply(.5, convert<svec>(vAdd(
			cheb::power(-1, range),1)));
		if (other(switch_[0]) > this->operator()(switch_[0]))
			keys = vSub(1.,keys);
		std::vector<IntervalFunction> funs;
		auto ints = switch_.intervals();
		for (int32_t i = 0; i < ints.size(); ++i) {
			auto subdom = newdom.restricted(Domain{ ints[i].a, ints[i].b });
			Function subfun;
			if (keys[i])
				subfun = restricted(subdom);
			else
				subfun = other.restricted(subdom);
			funs.insert(std::end(funs), std::begin(subfun.funs), std::end(subfun.funs));
		}
		return Function(funs);
	}
	Function Function::maximum(const scalar& c) const {
		return maximum(Function(c, domain()));
	}
	Function Function::minimum(const Function& other) const {
		if (isempty() || other.isempty()) return Function();
		auto roots = (*this - other).roots();
		auto newdom = domain().united(other.domain()).merged((roots));
		auto interval = newdom.support();
		auto switch_ = Domain{ interval }.merged((roots));
		auto range = convert<lvec>(arange(0, (len)switch_.size() - 1));
		auto keys = vMultiply(.5,cheb::convert<svec>(vAdd(
			cheb::power(-1, range)
			,1)));
		if (other(switch_[0]) < this->operator()(switch_[0]))
			keys = vSub(1.,keys);
		std::vector<IntervalFunction> funs;
		auto ints = switch_.intervals();
		for (int32_t i = 0; i < ints.size(); ++i) {
			auto subdom = newdom.restricted(Domain{ ints[i].a, ints[i].b });
			Function subfun;
			if (keys[i])
				subfun = restricted(subdom);
			else
				subfun = other.restricted(subdom);
			funs.insert(std::end(funs), std::begin(subfun.funs), std::end(subfun.funs));
		}
		return Function(funs);
	}
	Function Function::minimum(const scalar& c) const {
		return minimum(Function(c, domain()));
	}
	Function pow(const Function& lhs, const Function& rhs) {
		if (lhs.isempty() || rhs.isempty())
			return Function();
		auto newdom = lhs.domain().united(rhs.domain());
		auto chbfn1 = lhs.breakWith(newdom);
		auto chbfn2 = rhs.breakWith(newdom);
		std::vector<IntervalFunction> newfuns;
		for (int32_t i = 0; i < chbfn1.funs.size(); ++i)
			newfuns.push_back(pow(chbfn1.funs[i], chbfn2.funs[i]).simplified());
		return Function(newfuns);
	}
	Function pow(const Function& lhs, const scalar& rhs) {
		if (lhs.isempty())
			return Function();
		auto chbfn1 = lhs;
		std::vector<IntervalFunction> newfuns;
		for (int32_t i = 0; i < chbfn1.funs.size(); ++i)
			newfuns.push_back(pow(chbfn1.funs[i], rhs));
		return Function(newfuns);
	}
	Function pow(const scalar& lhs, const Function& rhs) {
		if (rhs.isempty())
			return Function();
		auto chbfn1 = rhs;
		std::vector<IntervalFunction> newfuns;
		for (int32_t i = 0; i < chbfn1.funs.size(); ++i)
			newfuns.push_back(pow(lhs, chbfn1.funs[i]));
		return Function(newfuns);
	}
}