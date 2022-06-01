// #pragma once 
#include <cheb/cheb.h>
namespace cheb {
	IntervalFunction::IntervalFunction(vfunc f, Interval interval) : minterval(interval) {
		auto uifunc = [f, interval](svec s) {
			return f(interval(s));
		};
		_coeffs = (adaptive(uifunc));
	}
	IntervalFunction::IntervalFunction(vfunc f, len n, Interval interval) : minterval(interval) {
		auto uifunc = [f, interval](svec s) {
			auto si = interval(s);
			auto yi = f(si);
			return yi;
		};
		auto points = (chebpts2(n));
		auto values = uifunc(points);
		_coeffs = (vals2coeffs2((values)));
	}
	svec IntervalFunction::operator()(svec x, method how) const {
		auto y = minterval.mapFromInterval(x);
		if (how == method::clenshaw)
			return clenshaw(y, _coeffs);
		else {
			auto fk = values();
			auto xk = (chebpts2((len)fk.size()));
			auto vk = (barywts2((len)fk.size()));
			return (bary((y), (fk), (xk), (vk)));
		}
	}
	svec IntervalFunction::eval(svec x, method how) const {
		auto y = x;
		if (how == method::clenshaw)
			return (clenshaw((y), (_coeffs)));
		else {
			auto fk = values();
			auto xk = (chebpts2((len)fk.size()));
			auto vk = (barywts2((len)fk.size()));
			return (bary((y), (fk), (xk), (vk)));
		}
	}
	svec IntervalFunction::values() const {
		return (coeffs2vals2((_coeffs)));
	}
	scalar IntervalFunction::operator()(scalar x, method how) const {
		return this->operator()(svec{ x }, how)[0];
	}
	svec IntervalFunction::call_clenshaw(svec x) const {
		return this->operator()(x, method::clenshaw);
	}
	svec IntervalFunction::call_bary(svec x) const {
		return this->operator()(x, method::bary);
	}
	svec IntervalFunction::endvalues() const {
		return this->operator()(support());
	}
	Interval IntervalFunction::interval() const { return minterval; }
	svec IntervalFunction::support() const { return svec{ minterval.a,minterval.b }; }
	svec IntervalFunction::coeffs() const { return _coeffs; }
	len IntervalFunction::size() const { return (len)_coeffs.size(); }
	bool IntervalFunction::isconst() const { return _coeffs.size() == 1; }
	bool IntervalFunction::isempty() const { return _coeffs.size() == 0; }
	scalar IntervalFunction::vscale() const {
		if (isempty()) return 0.0;
		auto vals = _coeffs;
		if (!isconst())
			vals = values();
		auto vscale = 0.0;
		for (auto v : vals)
			vscale = std::max(vscale, std::abs(v));
		return vscale;
	}
	IntervalFunction IntervalFunction::prolonged(int32_t n) const {
		auto m = (len)_coeffs.size();
		auto ak = _coeffs;
		if ((n - m) < 0)
			return IntervalFunction(subVec(ak, n), interval());
		else if ((n - m) > 0) {
			auto vk = (ak);
			for (int32_t i = 0; i < (n - m); ++i)
				vk.push_back(0.0);
			return IntervalFunction((vk), interval());
		}
		return IntervalFunction(_coeffs, interval());
	}

	IntervalFunction IntervalFunction::restricted(Interval subinterval) const {
		if (!minterval.contains(subinterval))
			throw std::invalid_argument("Not SubInterval");
		if (minterval == subinterval)
			return *this;
		return IntervalFunction(
			[*this](svec s) {return this->operator()(s); }, size(), subinterval);
	}
	svec IntervalFunction::roots(bool polish, int32_t splitDegree)const {
		auto roots = rootsunit(_coeffs, 1e2 * cheb::eps, splitDegree);
		if (!polish) return minterval(roots);
		auto nroots = (newtonroots(*this, (roots)));
		return minterval(nroots);
	}
	IntervalFunction IntervalFunction::antiDerivative() const {
		if (size() == 0) return *this;
		auto	n = size();
		auto	ak = concatenate(_coeffs, { 0, 0 });
		auto	bk = svec(n + 1, 0.);
		auto rk = arange(2, n + 1);
		if (rk.size() > 0) {
			for (int32_t i = 0; i < n - 1; ++i)
				bk[2 + i] = .5 * (ak[1 + i] - ak[3 + i]) / rk[i];
		}
		bk[1] = ak[0] - .5 * ak[2];
		auto vk = svec(n, 1.);
		for (int32_t i = 1; i < n; i += 2)
			vk[i] = -1;
		bk[0] = cheb::sum(vMultiply(vk, svec(bk.begin() + 1, bk.end())));
		auto cum = IntervalFunction(bk, minterval);
		auto [a, b] = minterval;
		return (.5 * (b - a)) * cum;
	}
	IntervalFunction IntervalFunction::derivative()const {
		if (size() == 0) return *this;
		if (isconst()) return IntervalFunction(0., minterval);
		auto n = size();
		auto ak = _coeffs;
		auto zk = svec(n - 1, 0.);
		auto wk = (vMultiply(2, arange(1, n)));
		auto vk = vMultiply(wk, svec(ak.begin() + 1, ak.end()));
		if (zk.size() > 0) {
			zk[n - 2] = vk[n - 2];
		}
		if (zk.size() > 1) {
			zk[n - 3] = vk[n - 3];
		}
		for (int32_t i = n - 4; i >= 0; i -= 2)
			zk[i] = vk[i] + zk[i + 2];
		for (int32_t i = n - 5; i >= 0; i -= 2)
			zk[i] = vk[i] + zk[i + 2];
		zk[0] = .5 * zk[0];
		IntervalFunction diff(zk, minterval);
		auto [a, b] = minterval;
		return (2. / (b - a)) * diff;
	}
	scalar IntervalFunction::defIntegral()const {
		if (_coeffs.size() == 0) return 0.;
		if (isconst())
			return 2. * this->operator()(0.);
		auto ak = _coeffs;
		for (int32_t i = 1; i < ak.size(); i += 2)
			ak[i] = 0;
		auto kk = arange(2, (len)ak.size());
		auto ii = concatenate(svec{ 2,0 }, vDivide(2., convert<svec>(vSub(1, vMultiply(kk, kk)))));
		auto sum = cheb::sum(vMultiply(ak, ii));

		auto [a, b] = minterval;
		return (.5 * (b - a)) * sum;
	}
	IntervalFunction IntervalFunction::simplified() const {
		auto cfs = _coeffs;
		auto n = standard_chop((cfs));
		return IntervalFunction(subVec(cfs, n), minterval);
	}
	bool IntervalFunction::operator==(const IntervalFunction& rhs) {
		return (minterval == rhs.minterval) && (_coeffs == rhs._coeffs);
	}
	IntervalFunction pow(const IntervalFunction& lhs, const IntervalFunction& rhs) {
		if (lhs.interval() != rhs.interval()) throw std::invalid_argument("interval mismatch");
		if (lhs.isempty() || rhs.isempty()) return IntervalFunction();
		return IntervalFunction([lhs, rhs](scalar x) {
			return std::pow(lhs(x), rhs(x));
			}, lhs.interval());
	}
	IntervalFunction pow(const IntervalFunction& lhs, const scalar& rhs) {
		if (lhs.isempty()) return IntervalFunction();
		return IntervalFunction([lhs, rhs](scalar x) {
			return std::pow(lhs(x), rhs);
			}, lhs.interval());
	}
	IntervalFunction pow(const scalar& lhs, const IntervalFunction& rhs) {
		return IntervalFunction(
			[lhs, rhs](svec x) { return cheb::power(lhs, rhs(x)); }
		, rhs.interval());
	}
}
using namespace cheb;

IntervalFunction operator+(const IntervalFunction& lhs, const IntervalFunction& rhs) {
	if (lhs.interval() != rhs.interval()) throw std::invalid_argument("interval mismatch"); \
		auto g = lhs;
	auto f = rhs;
	if (lhs.isempty())
		return lhs;
	if (f.isempty())
		return f;
	auto n = lhs.size();
	auto m = rhs.size();
	if (n < m)
		g = g.prolonged(m);
	else if (m < n)
		f = f.prolonged(n);
	svec coeffs(std::max(n, m));
	for (int32_t i = 0; i < std::max(n, m); ++i)
		coeffs[i] = g._coeffs[i] + f._coeffs[i];
	auto tol = .2 * eps * std::max(f.vscale(), g.vscale());
	if (all(coeffs, [tol](scalar v) {return abs(v) < tol; }))
		return IntervalFunction(0., lhs.interval());
	return IntervalFunction(coeffs, lhs.interval());
}
IntervalFunction operator+(const IntervalFunction& lhs, const scalar& rhs) {
	if (lhs.isempty())
		return lhs;
	auto cfs = lhs._coeffs;
	cfs[0] += rhs;
	return IntervalFunction(cfs, lhs.interval());
}
IntervalFunction operator+(const scalar& rhs, const IntervalFunction& lhs) {
	return lhs + rhs;
}
IntervalFunction operator-(const IntervalFunction& lhs, const IntervalFunction& rhs) {
	if (lhs.interval() != rhs.interval()) throw std::invalid_argument("interval mismatch"); \
		return lhs + (-rhs);
}
IntervalFunction operator-(const IntervalFunction& lhs, const scalar& rhs) {
	if (lhs.isempty())
		return lhs;
	return lhs + (-rhs);
}
IntervalFunction operator-(const scalar& lhs, const IntervalFunction& rhs) {
	return lhs + (-rhs);
}
IntervalFunction operator*(const IntervalFunction& lhs, const IntervalFunction& rhs) {
	if (lhs.interval() != rhs.interval()) throw std::invalid_argument("interval mismatch"); \
		if (lhs.isempty())
			return lhs;
	if (rhs.isempty())
		return rhs;
	auto g = lhs;
	auto f = rhs;
	auto n = f.size() + g.size() - 1;
	f = f.prolonged(n);
	g = g.prolonged(n);
	auto cfs = coeffmult(f.coeffs(), g.coeffs());
	return IntervalFunction(cfs, lhs.interval());
}
IntervalFunction operator*(const IntervalFunction& lhs, const scalar& rhs) {
	if (lhs.isempty())
		return lhs;
	auto cfs = evaluate(lhs._coeffs, [rhs](auto val) {return val * rhs; });
	return IntervalFunction(cfs, lhs.interval());

}
IntervalFunction operator*(const scalar& lhs, const IntervalFunction& rhs) {
	return rhs * lhs;
}
IntervalFunction operator/(const IntervalFunction& lhs, const IntervalFunction& rhs) {
	if (lhs.interval() != rhs.interval()) throw std::invalid_argument("interval mismatch"); \
		if (lhs.isempty())
			return lhs;
	if (rhs.isempty())
		return rhs;
	return IntervalFunction([lhs, rhs](scalar x) {return lhs(x) / rhs(x); }, lhs.interval());
}
IntervalFunction operator/(const IntervalFunction& lhs, const scalar& rhs) {
	return lhs * (1.0 / rhs);
}
IntervalFunction operator/(const scalar& lhs, const IntervalFunction& rhs) {
	if (rhs.isempty())
		return rhs;
	return IntervalFunction([lhs, rhs](scalar x) {return lhs / rhs(x); }, rhs.interval());
}
