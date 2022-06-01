#include <cheb/cheb.h>
#include <cheb/vectorutil.h>
#include <random>
#include <sstream>

namespace cheb {
	bool all(const bvec& cfs) {
		bool init = true;
		for (const auto& val : cfs) {
			init = init && val;
		}
		return init;
	}
	bool all(const svec& cfs, std::function<bool(scalar)> predicate) {
		return all(evaluate(cfs, predicate));
	}
	bool all(const cvec& cfs, std::function<bool(complex)> predicate) {
		return all(evaluate(cfs, predicate));
	}
	bool any(const bvec& cfs) {
		bool init = false;
		for (const auto& val : cfs) {
			init = init || val;
		}
		return init;
	}
	bool any(const svec& cfs, std::function<bool(scalar)> func) {
		return any(evaluate(cfs, func));
	}
	bool any(const cvec& cfs, std::function<bool(complex)> func) {
		return all(evaluate(cfs, func));
	}
	scalar infnorm(const svec& vals) {
		scalar m = 0.0;
		for (auto v : vals)
			m = std::max(std::abs(v), m);
		return m;
	}
	scalar infnorm(const svec& a, const svec& b) {
		return infnorm(vSub(a, b));
	}
	svec randn(int32_t n) {
		static std::random_device rd;
		static std::default_random_engine gen(rd());
		static std::uniform_real_distribution<scalar> dis(0.0, 1.0);
		svec randData(n);
		for (int32_t i = 0; i < n; ++i)
			randData[i] = (dis(gen));
		return randData;
	}
	svec linspace(scalar min, scalar max, int32_t n) {
		svec data(n);
		for (int32_t i = 0; i < n; ++i)
			data[i] = (min + (max - min) / ((scalar)n - 1.0) * (scalar)i);
		return data;
	}
	svec linspaceV(scalar min, scalar max, int32_t n) {
		svec data(n);
		for (int32_t i = 0; i < n; ++i)
			data[i] = (min + (max - min) / ((scalar)n - 1.0) * (scalar)i);
		return data;
	}
	lvec arange(len min, len max) {
		lvec valA(max - min,0);
		for (int32_t i = 0; i < max - min; ++i)
			valA[i] = min + i;
		return valA;
	}
	vfunc funcToVfunc(func fn) {
		return [fn](svec x) -> svec {return evaluate(x, fn); };
	}
}
