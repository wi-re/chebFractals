
//#define _CRT_SECURE_NO_WARNINGS
#include "SPH.h"
#include "2DMath.h"
#include "config.h"
#include <algorithm>
#include <array>
#include <boost/range/combine.hpp>
#include <iostream>
#include <numeric>
#include <chrono>
#include <sstream>
#include <atomic>


#include <tools/timer.h>
#include <cfloat>
#include <filesystem>

using clk = std::chrono::high_resolution_clock;
scalar toMs(clk::duration dur) {
  return static_cast<scalar>(std::chrono::duration_cast<std::chrono::microseconds>(dur).count()) / scalar(1000.0);
}
std::vector<std::tuple<cheb::complex, cheb::complex, cheb::complex>> trace;
std::filesystem::path expand(std::filesystem::path in) {
	namespace fs = std::filesystem;
#ifndef _WIN32
	if (in.string().size() < 1) return in;

	const char * home = getenv("HOME");
	if (home == NULL) {
		std::cerr << "error: HOME variable not set." << std::endl;
		throw std::invalid_argument("error: HOME environment variable not set.");
	}

	std::string s = in.string();
	if (s[0] == '~') {
		s = std::string(home) + s.substr(1, s.size() - 1);
		return fs::path(s);
	}
	else {
		return in;
	}
#else
	if (in.string().size() < 1) return in;

	const char * home = getenv("USERPROFILE");
	if (home == NULL) {
		std::cerr << "error: USERPROFILE variable not set." << std::endl;
		throw std::invalid_argument("error: USERPROFILE environment variable not set.");
	}

	std::string s = in.string();
	if (s[0] == '~') {
		s = std::string(home) + s.substr(1, s.size() - 1);
		return fs::path(s);
	}
	else {
		return in;
	}
#endif 

}
std::filesystem::path resolveFile(std::string fileName, std::vector<std::string> search_paths) {
	namespace fs = std::filesystem;
	auto& pm = ParameterManager::instance();
	fs::path working_dir = pm.get<std::string>("internal.working_directory");
	fs::path binary_dir = pm.get<std::string>("internal.binary_directory");
	fs::path source_dir = pm.get<std::string>("internal.source_directory");
	fs::path build_dir = pm.get<std::string>("internal.build_directory");
	fs::path expanded = expand(fs::path(fileName));

	fs::path base_path = "";
        if (fs::exists(expand(fs::path(fileName))))
          return expand(fs::path(fileName));
	for (const auto& path : search_paths){
		auto p = expand(fs::path(path));
		if (fs::exists(p / fileName))
			return p.string() + std::string("/") + fileName;
}

	if (fs::exists(fileName)) return fs::path(fileName);
        if (fs::exists(expanded))
          return expanded;

        for (const auto &pathi : search_paths) {
                        auto path = expand(fs::path(pathi));
          if (fs::exists(working_dir / path / fileName))
            return (working_dir / path / fileName).string();
          if (fs::exists(binary_dir / path / fileName))
            return (binary_dir / path / fileName).string();
          if (fs::exists(source_dir / path / fileName))
            return (source_dir / path / fileName).string();
          if (fs::exists(build_dir / path / fileName))
            return (build_dir / path / fileName).string();
        }

	if (fs::exists(working_dir / fileName))
          return (working_dir / fileName);
        if (fs::exists(binary_dir / fileName))
          return (binary_dir / fileName);
        if (fs::exists(source_dir / fileName))
          return (source_dir / fileName);
        if (fs::exists(build_dir / fileName))
          return (build_dir / fileName);

	std::stringstream sstream;
	sstream << "File '" << fileName << "' could not be found in any provided search path" << std::endl;
	std::cerr << sstream.str();
	throw std::runtime_error(sstream.str().c_str());
}


cheb::complex clenshaw(cheb::complex x, const cheb::svec& ak) {
	if (ak.size() == 0)
		return 0.0;
	if (ak.size() == 1)
		return ak[0];
	//if (any(ak, [](auto s) {return std::isnan(s); }))
	//    return std::nan("") * svec(1., xx.size());
	cheb::complex bk1 = 0.0;
	cheb::complex bk2 = 0.0;
	x *= 2.0;

	for (int32_t k = (int32_t)ak.size() - 1; k > 1; k -= 2) {
		bk2 = ak[k] + x * bk1 - bk2;
		bk1 = ak[k - 1] + x * bk2 - bk1;
	}
	if ((ak.size() - 1) % 2 == 1) {
		auto temp = bk1;
		bk1 = ak[1] + x * bk1 - bk2;
		bk2 = temp;
	}
	cheb::complex res = 0.0;
	res = ak[0] + 0.5 * x * bk1 - bk2;

	return res;
}

	struct complexState{
		FunctionState r, i;
	};

	complexState clenshawDeriv(cheb::complex val, const cheb::svec& ak) {
		if (ak.size() == 0)
			return complexState{};
		if (ak.size() == 1)
			return complexState{.r = FunctionState{.f = ak[0]}};

		std::vector<complexState> b(ak.size() + 2);
		auto x = val.real();
		auto y = val.imag();


		for (int32_t k = (int32_t)ak.size() - 1; k >= 0; k -= 1) {	
			auto s = k == 0 ? 1. : 2.;


			b[k].r.f        = ak[k] - b[k+2].r.f        + s * (x * b[k+1].r.f                                            - y * b[k+1].i.f);

			b[k].r.J.dfdx   =       - b[k+2].r.J.dfdx   + s * (x * b[k+1].r.J.dfdx                     + b[k+1].r.f      - y * b[k+1].i.J.dfdx);

			b[k].r.H.d2fdx2 =       - b[k+2].r.H.d2fdx2 + s * (x * b[k+1].r.H.d2fdx2 + b[k+1].r.J.dfdx + b[k+1].r.J.dfdx - y * b[k+1].i.H.d2fdx2);
			b[k].r.H.d2fdxy =       - b[k+2].r.H.d2fdxy + s * (x * b[k+1].r.H.d2fdxy + b[k+1].r.J.dfdy                   - y * b[k+1].i.H.d2fdxy - b[k+1].i.J.dfdx);

			b[k].r.J.dfdy   =       - b[k+2].r.J.dfdy   + s * (x * b[k+1].r.J.dfdy                                       - y * b[k+1].i.J.dfdy - b[k+1].i.f);

			b[k].r.H.d2fdyx =       - b[k+2].r.H.d2fdyx + s * (x * b[k+1].r.H.d2fdyx + b[k+1].r.J.dfdy                   - y * b[k+1].i.H.d2fdyx - b[k+1].i.J.dfdx);
			b[k].r.H.d2fdy2 =       - b[k+2].r.H.d2fdy2 + s * (x * b[k+1].r.H.d2fdy2                                     - y * b[k+1].i.H.d2fdy2 - b[k+1].i.J.dfdy - b[k+1].i.J.dfdy);

			b[k].i.f        =       - b[k+2].i.f        + s * (x * b[k+1].i.f                                            + y * b[k+1].r.f);

			b[k].i.J.dfdx   =       - b[k+2].i.J.dfdx   + s * (x * b[k+1].i.J.dfdx                     + b[k+1].i.f      + y * b[k+1].r.J.dfdx);

			b[k].i.H.d2fdx2 =       - b[k+2].i.H.d2fdx2 + s * (x * b[k+1].i.H.d2fdx2 + b[k+1].i.J.dfdx + b[k+1].i.J.dfdx + y * b[k+1].r.H.d2fdx2);
			b[k].i.H.d2fdxy =       - b[k+2].i.H.d2fdxy + s * (x * b[k+1].i.H.d2fdxy + b[k+1].i.J.dfdy                   + y * b[k+1].r.H.d2fdxy + b[k+1].r.J.dfdx);

			b[k].i.J.dfdy   =       - b[k+2].i.J.dfdy   + s * (x * b[k+1].i.J.dfdy                                       + y * b[k+1].r.J.dfdy   + b[k+1].r.f);

			b[k].i.H.d2fdyx =       - b[k+2].i.H.d2fdyx + s * (x * b[k+1].i.H.d2fdyx + b[k+1].i.J.dfdy                   + y * b[k+1].r.H.d2fdyx + b[k+1].r.J.dfdx);
			b[k].i.H.d2fdy2 =       - b[k+2].i.H.d2fdy2 + s * (x * b[k+1].i.H.d2fdy2                                     + y * b[k+1].r.H.d2fdy2 + b[k+1].r.J.dfdy + b[k+1].r.J.dfdy);
		}

		return b[0];
	}


cheb::Function globalFunction;
cheb::Function globalFunctionFirstDerivative;
cheb::Function globalFunctionSecondDerivative;

cheb::complex evalFunction(cheb::complex location){
    return clenshaw(location, globalFunction.funs[0].coeffs());
}
std::mutex m;
cheb::complex evalDerivative(cheb::complex location){
	auto [fr, fi] = clenshawDeriv(location, globalFunction.funs[0].coeffs());

	cheb::complex f(fr.f, fi.f);
	Jacobian J(fr.J.dfdx, fi.J.dfdx);
	Hessian H(fr.H.d2fdx2, fr.H.d2fdxy, fi.H.d2fdx2, fi.H.d2fdxy);

	// {
	// 	std::lock_guard lg(m);
	// 	// auto fc = clenshaw(location, globalFunction.funs[0].coeffs());
	// 	// auto hx = 1e-8;
	// 	// auto x = clenshaw(location + cheb::complex(0., 0.), globalFunction.funs[0].coeffs());
	// 	// auto xr = clenshaw(location + cheb::complex(hx, 0.), globalFunction.funs[0].coeffs());
	// 	// auto xi = clenshaw(location + cheb::complex(0., hx), globalFunction.funs[0].coeffs());

	// 	// std::cout << location.real() << " : " << location.imag() << " -> " << fx.real() << " : " << fx.imag() << " - " << fc.real() << " - " << fc.imag() << std::endl;
	// 	// std::cout << "Clenshaw: " << c.real() << " : " << c.imag() << std::endl;
	// 	// std::cout << "Jacobian: " << Jx.first.real() << " : " << Jx.second.real() << " \\ " << Jx.first.imag() << " : " << Jx.second.imag() << std::endl;
	// 	// std::cout << "Finite:   " << (xr - x).real() / hx << " : " << (xr - x).imag() / hx << " \\ " << (xi - x).real() / hx << " : " << (xi - x).imag() / hx << std::endl;
	// 	// std::cout << std::endl;
	// 	// std::cout << c.real() / Jx.first.real() << " : " << c.imag() / Jx.second.real() << std::endl;
	// 	// std::cout << c.real() / Jx.first.imag() << " : " << c.imag() / Jx.second.imag() << std::endl;
	// 	std::cout << "Evaluated at x = " << location.real() << " + " << location.imag() << "i" << std::endl;
	// 	std::cout << "Analytical:\n";

	// 	std::cout << "\t f(x) = " << fr.f << " + " << fi.f << "i\n";
	// 	std::cout << "\t J(x) = [ " << fr.J.dfdx << " " << fr.J.dfdy << " ] \\\\ [ " << fi.J.dfdx << " " << fi.J.dfdy << "]\n";
	// 	std::cout << "\t H(x)r= [ " << fr.H.d2fdx2 << " " << fr.H.d2fdxy << " ] \\\\ [ " << fr.H.d2fdyx << " " << fr.H.d2fdy2 << "]\n";
	// 	std::cout << "\t H(x)i= [ " << fi.H.d2fdx2 << " " << fi.H.d2fdxy << " ] \\\\ [ " << fi.H.d2fdyx << " " << fi.H.d2fdy2 << "]\n";
		
	// 	auto hx = 1e-8;
	// 	auto x   = clenshaw(location + cheb::complex(0., 0.), globalFunction.funs[0].coeffs());
	// 	auto xrp = clenshaw(location + cheb::complex(hx, 0.), globalFunction.funs[0].coeffs());
	// 	auto xip = clenshaw(location + cheb::complex(0., hx), globalFunction.funs[0].coeffs());
	// 	auto xrm = clenshaw(location - cheb::complex(hx, 0.), globalFunction.funs[0].coeffs());
	// 	auto xim = clenshaw(location - cheb::complex(0., hx), globalFunction.funs[0].coeffs());

	// 	std::cout << "\nNumerical:\n";
	// 	std::cout << "Finite:   " << (xrp - x).real() / hx << " : " << (xip - x).real() / hx << " \\ " << (xrp - x).imag() / hx << " : " << (xip - x).imag() / hx << std::endl;


	// 	std::cout << "\nClenshaw:\n";
	// 	auto dc = clenshaw(location, globalFunctionFirstDerivative.funs[0].coeffs());
	// 	auto d2c = clenshaw(location, globalFunctionSecondDerivative.funs[0].coeffs());
	// 	std::cout << "\t f'(x) = "  << dc.real() << " + " << dc.imag() << "i\n";
	// 	std::cout << "\t f''(x) = " << d2c.real() << " + " << d2c.imag() << "i\n";
	// 	{
			
	// 		auto hx = 1e-8;
	// 		auto x   = clenshaw(location + cheb::complex(0., 0.), globalFunctionFirstDerivative.funs[0].coeffs());
	// 		auto xrp = clenshaw(location + cheb::complex(hx, 0.), globalFunctionFirstDerivative.funs[0].coeffs());
	// 		auto xip = clenshaw(location + cheb::complex(0., hx), globalFunctionFirstDerivative.funs[0].coeffs());
	// 		auto xrm = clenshaw(location - cheb::complex(hx, 0.), globalFunctionFirstDerivative.funs[0].coeffs());
	// 		auto xim = clenshaw(location - cheb::complex(0., hx), globalFunctionFirstDerivative.funs[0].coeffs());

	// 		std::cout << "\nNumerical:\n";
	// 		std::cout << "Finite Hessian:   " << (xrp - x).real() / hx << " : " << (xip - x).real() / hx << " \\ " << (xrp - x).imag() / hx << " : " << (xip - x).imag() / hx << std::endl;
	// 	}


	// 	std::cout << std::endl << std::endl;
	
	
	// }
	return cheb::complex(fr.J.dfdx, fi.J.dfdx);
	//return cheb::complex(Jx.second.imag(), -Jx.second.real());

    return clenshaw(location, globalFunctionFirstDerivative.funs[0].coeffs());
}
cheb::complex evalSecondDerivative(cheb::complex location){
    return clenshaw(location, globalFunctionSecondDerivative.funs[0].coeffs());
}


std::tuple<cheb::complex, cheb::complex, cheb::complex> evalPolynomial(cheb::complex location){
	auto [fr, fi] = clenshawDeriv(location, globalFunction.funs[0].coeffs());

	cheb::complex f(fr.f, fi.f);
	Jacobian J(fr.J.dfdx, fi.J.dfdx);
	Hessian H(fr.H.d2fdx2, fr.H.d2fdxy, fi.H.d2fdx2, fi.H.d2fdxy);
	return std::make_tuple(f, cheb::complex(J.dfdx, J.dfdy), cheb::complex(H.d2fdx2,H.d2fdyx));
}

FunctionState evalSquarePolynomial(cheb::complex location){
	auto [fr, fi] = clenshawDeriv(location, globalFunction.funs[0].coeffs());

	auto f = fr.f * fr.f + fi.f * fi.f;

	auto dfdx = 2. * fr.f * fr.J.dfdx + 2. * fi.f * fi.J.dfdx;
	auto dfdy = 2. * fr.f * fr.J.dfdy + 2. * fi.f * fi.J.dfdy;

	auto d2fdx2 = 2. * (fr.J.dfdx * fr.J.dfdx + fr.f * fr.H.d2fdx2 + fi.J.dfdx * fi.J.dfdx + fi.f * fi.H.d2fdx2);
	auto d2fdxy = 2. * (fr.J.dfdy * fr.J.dfdx + fr.f * fr.H.d2fdxy + fi.J.dfdy * fi.J.dfdx + fi.f * fi.H.d2fdxy);
	auto d2fdyx = 2. * (fr.J.dfdy * fr.J.dfdx + fr.f * fr.H.d2fdyx + fi.J.dfdy * fi.J.dfdx + fi.f * fi.H.d2fdyx);
	auto d2fdy2 = 2. * (fr.J.dfdy * fr.J.dfdy + fr.f * fr.H.d2fdy2 + fi.J.dfdy * fi.J.dfdy + fi.f * fi.H.d2fdy2);

	Jacobian J{dfdx,dfdy};
	Hessian H{d2fdx2, d2fdxy, d2fdyx, d2fdy2};

	return FunctionState{f,J,H};
}
std::vector<std::vector<double>> globalCoefficients;