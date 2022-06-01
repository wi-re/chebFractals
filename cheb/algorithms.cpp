#include <cheb/cheb.h>
#include <cheb/algorithms.h>
#include <armadillo>
#include <numbers>
#define FFTW_DLL 
#include <fftw3.h>
#include <array>
namespace cheb {
//#define USE_FFTW3
#ifdef USE_FFTW3
	std::mutex fftwlock;
	template<std::size_t N = 131072>
	struct fft_wrapper {
		fftw_complex* in, *out;
		std::map<std::size_t, fftw_plan> iplans, fplans;


		svec ifft(svec vals) {
			auto n = vals.size();
			if (!iplans.contains(n)) {
				std::lock_guard lg(fftwlock);
				iplans[n] = fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
			}
			
			for (int i = 0; i < n; i++) {
				in[i][0] = vals[i];
				in[i][1] = 0.0;
			}
			fftw_execute(iplans[n]);
			svec resv(vals.size());
			for (int i = 0; i < n; i++) {
				resv[i] = 1.0 / (double) n * out[i][0];
			}
			return resv;
		}
		svec fft(svec vals) {
			auto n = vals.size();
			if (!fplans.contains(n)) {
				std::lock_guard lg(fftwlock);
				fplans[n] = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
			}

			for (int i = 0; i < n; i++) {
				in[i][0] = vals[i];
				in[i][1] = 0.0;
			}
			fftw_execute(fplans[n]);
			svec resv(vals.size());
			for (int i = 0; i < n; i++) {
				resv[i] = out[i][0];
			}
			return resv;
		}
		~fft_wrapper() {
			std::lock_guard lg(fftwlock);
			for (auto& p : iplans)
				fftw_destroy_plan(p.second);
			for (auto& p : fplans)
				fftw_destroy_plan(p.second);
			delete[] in;
			delete[] out;
		}
		fft_wrapper() {
			in = new fftw_complex[N];
			out = new fftw_complex[N];
		}
	};
	thread_local fft_wrapper fftw;
#endif
	inline svec fft(svec vals) {
		if (vals.size() <= 1)
			return vals;

#ifdef USE_FFTW3
		auto n = vals.size();
		if ((n & (n - 1)) == 0) {
			//fft_wrapper fftw;
			return fftw.fft(vals);
		}
		else 
#endif
		{
			arma::vec v = arma::zeros(vals.size());
			for (int32_t i = 0; i < vals.size(); ++i)
				v[i] = vals[i];
			arma::cx_vec res = arma::fft(v);
			svec resv(res.size(),0.);
			for (int32_t i = 0; i < res.size(); ++i)
				resv[i] = res[i].real();
			return resv;
		}

	}

	inline svec ifft(svec vals) {
		if (vals.size() <= 1)
			return vals;

#ifdef USE_FFTW3
		auto n = vals.size();
		if ((n & (n - 1)) == 0) {
			//fft_wrapper fftw;
			return fftw.ifft(vals);
		}
		else 
#endif
		{
			arma::cx_vec v(vals.size());
			for (int32_t i = 0; i < vals.size(); ++i)
				v[i] = vals[i];
			arma::cx_vec res = arma::ifft(v);
			svec resv(res.size(),0.);
			for (int32_t i = 0; i < res.size(); ++i)
				resv[i] = res[i].real();
			return resv;
		}
	}

	svec rootsunit(svec ak, scalar htol, int32_t splitDegree) {
		//return cheb::convert<svec>(rootsunit(cheb::convert<vscalar>(ak), htol));
		auto n = standard_chop((ak));
		ak = subVec(ak, n);

		// Recurse for large coefficient sets
		if (n > splitDegree) {
			auto chebpts = (chebpts2((len)ak.size()));
			auto lmap = Interval{ -1.,SPLITPOINT };
			auto rmap = Interval{ SPLITPOINT, 1. };
			auto lpts = lmap(chebpts);
			auto rpts = rmap(chebpts);
			auto lval = (clenshaw((lpts), (ak)));
			auto rval = (clenshaw((rpts), (ak)));
			auto lcfs = (vals2coeffs2((lval)));
			auto rcfs = (vals2coeffs2((rval)));
			auto lrts = rootsunit(lcfs, 2 * htol,splitDegree);
			auto rrts = rootsunit(rcfs, 2 * htol, splitDegree);
			auto ls = lmap((lrts));
			auto rs = rmap((rrts));
			return (concatenate(ls,rs));
		}
		// A degree of 1 or less polynomial has no roots.
		// Note that this implicitly means that f(x) = 0 has no roots.
		if (n <= 1)
			return svec();
		cvec rts;
		// If the number of coefficients is 2 we can directly find the root
		// as the function described here is a linear function 
		if (n == 2)
			rts = cvec{ -ak[0] / ak[1] };
		// Otherweise determine the roots as the real eigenvalues of the
		// Chebyshev Companion Matrix
		else if (n <= splitDegree) {
			arma::mat C = arma::zeros(n - 1, n - 1);
			C.diag(1) = arma::ones(n - 2) * .5;
			C.diag(-1) = arma::ones(n - 2) * .5;
			C(0, 1) = 1.;
			arma::mat D = arma::zeros(n - 1, n - 1);
			for (int32_t i = 0; i < n - 1; ++i) 
				D(n - 2, i) = ak[i];
			arma::mat E = C - .5 * 1. / (*(std::end(ak) - 1)) * D;

			arma::cx_vec eigval;
			//arma::cx_mat eigvec;
			arma::eig_gen(eigval, E);

			rts = cvec(eigval.size(),0.);
			for (int32_t i = 0; i < eigval.size(); ++i)
				rts[i] = eigval(i);
		}
		// First we need to remove all roots with a complex part larger than some threshold
		svec rrts;
		for (auto c : rts) {
			if (std::abs(c.imag()) < htol && (std::abs(c) <= 1. + 1. * htol))
				rrts.push_back(c.real());
		}
		//auto mask = std::abs(evaluate(rts, [](auto c) {return c.imag(); })) < htol;
		//auto rrts = evaluate(vcomplex(rts[mask]), [](auto c) {return c.real(); });
		// We then also remove any roots outside of the domain -1 1 (within some tolerance)
		// as these roots are not useful values and do not reflect actual roots.
		// Note that depending on the accuracy of the eigen solver, i.e. Eigen's eigen solver
		// this threshold needs to be fairly large to not miss any roots
		//svec rrts2 = rrts[std::abs(rrts) <= 1. + 1. * htol];
		// Sort the roots from left to right for convienence
		rrts = sort(rrts);
		// Ensure that the first and last root are in the correct interval to help the 
		// Newton polishing later
		if (rrts.size() >= 1) {
			rrts[0] = std::min(std::max(rrts[0], -1.), 1.);
			*(std::end(rrts) - 1) = std::min(*(std::end(rrts) - 1), 1.);
		}
		return rrts;
	}
	svec bary(svec xx, svec fk, svec xk, svec vk) {
		if (xx.size() == 0 || fk.size() == 0)
			return svec{};
		if (fk.size() == 1)
			return vMultiply(fk[0], svec(xx.size(),1.));
		for(auto v : fk)
			if(std::isnan(v))
				return vMultiply(std::nan(""), svec(xx.size(),1.));
		cheb::svec out;
		// If the number of evaluation points is less than the number of nodes
		// we iterate over the evaluation points
		if (xx.size() < 4 * xk.size()) {
			out = cheb::svec(xx.size()),0.0;
			for (int32_t i = 0; i < xx.size(); ++i) {
				auto tt = vDivide(vk, vSub(xx[i], xk));
				out[i] = dot(tt, fk) / cheb::sum(tt);
			}
		}
		// else iterate over the nodes
		else {
			auto numer = cheb::svec(xx.size(),0.);
			auto denom = cheb::svec(xx.size(),0.);
			for (int32_t j = 0; j < xk.size(); ++j) {
				//auto temp = vk[j] / (xx - xk[j]);
				auto temp = vDivide(vk[j], vSub(xx, xk[j]));
				numer = vAdd(numer, vMultiply(temp,fk[j]));
				denom = vAdd(denom, temp);
			}
			out = vDivide(numer, denom);
		}
		// In this process NaNs can arise, especially due to 
		// auto tt = vk / (xx[i] - xk);
		// These NaN's are always on function evaluation points that
		// coincide with the nodes of the polynomial. These NaNs can be resolved,
		// however, by setting the function values at these points to the function
		// values corresponding to the function values at the nodes (which can be
		// determined directly through the coefficients of the polynomial)
		for (int32_t k = 0; k < out.size(); ++k)
			if (std::isnan(out[k])) {
				auto idx = -1;
				for (int32_t i = 0; i < xk.size(); ++i)
					if (xx[k] == xk[i])
						idx = i;
				if (idx != -1)
					out[k] = fk[idx];
			}
		return out;
	}

	cheb::svec clenshaw(cheb::svec xx, const cheb::svec& ak) {
		if (xx.size() == 0)
			return std::vector<double>{};
		if (ak.size() == 0)
			return std::vector<double>(xx.size(), 0.0);
		if (ak.size() == 1)
			return std::vector<double>(xx.size(), ak[0]);
		//if (any(ak, [](auto s) {return std::isnan(s); }))
		//    return std::nan("") * svec(1., xx.size());
		auto bk1 = std::vector < double>(xx.size(), 0.0);
		auto bk2 = std::vector < double>(xx.size(), 0.0);
		for (auto& x : xx)
			x *= 2.0;
		for (int32_t k = (int32_t)ak.size() - 1; k > 1; k -= 2) {
			for (int32_t i = 0; i < xx.size(); ++i) {
				bk2[i] = ak[k] + xx[i] * bk1[i] - bk2[i];
				bk1[i] = ak[k - 1] + xx[i] * bk2[i] - bk1[i];
			}
		}
		if ((ak.size() - 1) % 2 == 1) {
			for (int32_t i = 0; i < xx.size(); ++i) {
				auto temp = bk1[i];
				bk1[i] = ak[1] + xx[i] * bk1[i] - bk2[i];
				bk2[i] = temp;
			}
		}
		std::vector<double> res(xx.size(), 0.0);
		for (int32_t i = 0; i < xx.size(); ++i)
			res[i] = ak[0] + 0.5 * xx[i] * bk1[i] - bk2[i];
		return res;
	}


	len round(scalar x)	{
		return (int)(x + 0.5);
	}
	len standard_chop(svec coeffs, scalar tol) {
		// This process only works if there are at least 17 coefficients
		auto n = (len)coeffs.size();
		auto cutoff = n;
		if (n < 17)
			return cutoff;
		// First convert the provided coefficients to a new vector envelope
		// that's monotically non-increasing (found using a scan of the flipped abs coefficients)
		auto b = coeffs;
		for (int32_t i = 0; i < n; ++i)
			b[n - 1 - i] = std::abs(coeffs[i]);
		//auto b = flip(std::abs(coeffs));
		auto m = flip(scan(b, [](auto l, auto r) {return std::max(l, r); }));
		// If m[0] == 0. we cannot normalize, return 1 instead as fallback
		if (m[0] == 0.) return 1;
		// Normalize to begin at 0
		auto envelope = vDivide(m,m[0]);
		// Find the first plateau point, i.e. the point at which the coefficients are non-increasing
		auto plateauPoint = 0;
		int32_t j2 = 0;
		for (int32_t j = 1; j < n; ++j) {
			j2 = round(1.25 * (double)j + 5.);
			if (j2 > n - 1) // There is no plateau so all coefficients are required
				return cutoff;
			auto e1 = envelope[j];
			auto e2 = envelope[j2];
			auto r = 3.0 * (1. - std::log(e1) / std::log(tol));
			//if (e1 == 0. || e2 == 0.)
			//	continue;
			bool plateau = (e1 == 0.) || (e2 / e1 > r);
			if (plateau) { // found the first plateau point
				//std::clog << e1 << " : " << e2 << "  -> " << e2 / e1 << " : " << r << " | " << tol << " . "  << std::log(e1) << " | " << std::log(tol) << std::endl;
				plateauPoint = j;
				break;
			}
		}
		//debugPrintVec(envelope);
		// The cutoff is the fixed to a point where the envelope, plus a linear function
		// included to bias the results towards the left end, is minimal.
		if (envelope[plateauPoint] == 0.)
			cutoff = plateauPoint;
		else {
			auto j3 = 0;
			for (int32_t i = 0; i < envelope.size(); ++i) {
				j3 += envelope[i] >= std::pow(tol, 7. / 6.) ? 1 : 0;
			}//evaluate(envelope >= std::pow(tol, 7. / 6.), [](auto v) {return v ? 1 : 0; }).sum();
			if (j3 < j2) {
				j2 = j3 + 1;
				envelope[j2] = std::pow(tol, 7. / 6.);
			}
			envelope.resize(j2);
			auto cc = vAdd(cheb::log10(envelope), (linspaceV(0, (-1. / 3.) * std::log10(tol),j2)));
			
			//cc += linspace(0, (-1. / 3.) * std::log10(tol), j2);
			auto d = argmin(cc);
			cutoff = d;
		}
		return std::min(cutoff, n - 1);
	}
	svec adaptive(vfunc fun, uint32_t maxpow2) {
		svec coeffs;
		for (int32_t k = 4; k < (len) maxpow2 + 1; ++k) {
			auto n = (len)std::pow(2, k) + 1;
			auto points = (chebpts2(n));
			auto values = fun(points);
			coeffs = vals2coeffs2((values));

			auto chplen = standard_chop(coeffs);

			//auto midpoints = svec(0., n - 1);
			//for (int32_t i = 0; i < n - 1; ++i)
			//	midpoints[i] = (points[i] + points[i + 1]) * .5;
			//auto cvals = clenshaw(midpoints, coeffs);
			//auto fvals = fun(midpoints);
			//auto error = infnorm(fvals - cvals);
			//auto tol = 1e3 * eps;
			//if (error <= tol) {
			//	return slice(coeffs, 0, standard_chop(coeffs));
			//}

			if (chplen < coeffs.size()) {
				coeffs.resize(chplen);
				return coeffs;
			}
				//return coeffs;
				//return slice(coeffs, 0, chplen);
			if (k == maxpow2) {
				throw(std::invalid_argument("Constructor did not converge"));
				std::cerr << "Constructor did not converge using " << n << " points" << std::endl;
				return coeffs;
			}
		}
		return coeffs;
	}
	svec newtonroots(cheb::IntervalFunction fun, svec rts, scalar tol, uint32_t maxiter) {
		auto [a, b] = fun.interval();
		if (rts.size() > 0) {
			auto dfun = fun.derivative() / (2. / (b - a));
			auto prv = vMultiply(std::numeric_limits<scalar>::infinity(), rts);
			auto count = 0;
			while (infnorm(vSub(rts, prv)) > tol && count++ <= (len)maxiter) {
				prv = rts;
				rts = vSub(rts, vDivide((fun.eval((rts))), (dfun.eval((rts)))));
			}
		}
		return rts;
	}

	svec coeffmult(svec fvc, svec gvc) {
		auto nf = (int32_t)fvc.size();
		auto ng = (int32_t)gvc.size();;
		svec Fc(1 + nf - 1 + nf - 1, 0.);
		Fc[0] = fvc[0] * 2.;
		for (int32_t i = 0; i < nf - 1; ++i) {
			Fc[1 + i] = fvc[1 + i];
			Fc[1 + i + nf - 1] = fvc[nf - 1 - i];
		}
		svec Gc(1 + ng - 1 + ng - 1, 0.);
		Gc[0] = gvc[0] * 2.;
		for (int32_t i = 0; i < ng - 1; ++i) {
			Gc[1 + i] = gvc[1 + i];
			Gc[1 + i + ng - 1] = gvc[ng - 1 - i];
		}

		auto ak = ifft(cheb::vMultiply(fft(Fc) ,fft(Gc)));

		svec result(ak.size() , 0.);
		result[0] = ak[0] * 0.25;
		auto ni = (int32_t)ak.size();
		for (int32_t i = 1; i < ni; ++i) {
			result[i] = (ak[i] + ak[ni-i]) * 0.25;
		}
		result.resize(nf);
		return result;
	}
	svec barywts2(len n) {
		if (n == 0)
			return svec{};
		if (n == 1)
			return svec{ 1. };
		svec wts(n,1);
		*(std::end(wts) - 1) = .5;
		for (int32_t i = n - 2; i >= 0; i -= 2)
			wts[i] = -1;
		wts[0] = .5 * wts[0];
		return wts;
	}
	svec chebpts2(len n) {
		if (n == 1)
			return svec{ 0. };
		auto nn = cheb::convert<svec>(flip((arange(0, n))));
		auto pts = cheb::cos(vMultiply(nn, std::numbers::pi_v<scalar> / ((scalar)n - 1.)));
		return pts;
	}
	svec vals2coeffs2(svec vals) {
		auto n = (len)vals.size();
		if (n <= 1) return vals;

		svec tmp(n + n - 2);
		for (int32_t i = 0; i < n; ++i) 
			tmp[i] = vals[n - 1 - i];

		for (int32_t i = 1; i < n-1; ++i) 
			tmp[n + i-1] = vals[i];
	
		//tmp[n - 1] = vals[0];

		// auto tmp = concatenate(flip(vals), slice(vals, 1, (len)vals.size() - 1));
		svec coeffs;
		coeffs = (ifft((tmp)));
		coeffs.resize(n);
		//coeffs = slice(coeffs, 0, n);
		for (int32_t i = 1; i < n - 1; ++i) {
			coeffs[i] = 2. * coeffs[i];
		}
		return (coeffs);
	}
	svec coeffs2vals2(svec coeffs) {
		auto n = (len)coeffs.size();
		if (n <= 1)
			return coeffs;
		for (int32_t i = 1; i < n - 1; ++i)
			coeffs[i] = .5 * coeffs[i];
		coeffs.resize(n + n - 1);


		for (int32_t i = 0; i < n-2; ++i)
			coeffs[n + i] = coeffs[n-2-i];

		//auto tmp = concatenate(coeffs, slice(coeffs, n - 2, 0, -1));
		svec values;
		values = (fft((coeffs)));
		values.resize(n);
		cheb::reverse(values);
		//values = slice(values, n - 1, -1, -1);
		return values;
	}

}