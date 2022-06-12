
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

	std::tuple<cheb::complex,std::pair<cheb::complex,cheb::complex>> clenshawDeriv(cheb::complex val, const cheb::svec& ak) {
		if (ak.size() == 0)
			return std::make_tuple(0.0,std::make_pair(0.0,0.0));
		if (ak.size() == 1)
			return std::make_tuple(ak[0], std::make_pair(0.0, 0.0));

		cheb::svec bkr(ak.size() + 2, 0.0);
		cheb::svec bki(ak.size() + 2, 0.0);
		cheb::svec bkr_x(ak.size() + 2, 0.0);
		cheb::svec bki_x(ak.size() + 2, 0.0);
		cheb::svec bkr_y(ak.size() + 2, 0.0);
		cheb::svec bki_y(ak.size() + 2, 0.0);

		auto x = val.real();
		auto y = val.imag();


		for (int32_t k = (int32_t)ak.size() - 1; k >= 1; k -= 1) {			
			bkr  [k] = ak[k] - bkr  [k + 2] + 2. * (x * bkr  [k + 1]              - y * bki  [k + 1]);
			bkr_x[k] =       - bkr_x[k + 2] + 2. * (x * bkr_x[k + 1] + bkr[k + 1] - y * bki_x[k + 1]);
			bkr_y[k] =       - bkr_y[k + 2] + 2. * (x * bkr_y[k + 1]              - y * bki_y[k + 1] - bki[k + 1]);

			bki  [k] =       - bki  [k + 2] + 2. * (x * bki  [k + 1]              + y * bkr  [k + 1]);
			bki_x[k] =       - bki_x[k + 2] + 2. * (x * bki_x[k + 1] + bki[k + 1] + y * bkr_x[k + 1]);
			bki_y[k] =       - bki_y[k + 2] + 2. * (x * bki_y[k + 1]              + y * bkr_y[k + 1] + bkr[k + 1]);
		}
		int32_t k = 0;
		bkr  [k] = ak[k] - bkr  [k + 2] + 1. * (x * bkr  [k + 1]              - y * bki  [k + 1]);
		bkr_x[k] =       - bkr_x[k + 2] + 1. * (x * bkr_x[k + 1] + bkr[k + 1] - y * bki_x[k + 1]);
		bkr_y[k] =       - bkr_y[k + 2] + 1. * (x * bkr_y[k + 1]              - y * bki_y[k + 1] - bki[k + 1]);

		bki  [k] =       -bki  [k + 2] + 1. * (x * bki  [k + 1]               + y * bkr  [k + 1]);
		bki_x[k] =       -bki_x[k + 2] + 1. * (x * bki_x[k + 1] + bki[k + 1]  + y * bkr_x[k + 1]);
		bki_y[k] =       -bki_y[k + 2] + 1. * (x * bki_y[k + 1]               + y * bkr_y[k + 1] + bkr[k + 1]);

		return std::make_tuple(
			cheb::complex(bkr[0],bki[0]), 
			std::make_pair(cheb::complex(bkr_x[0], bki_x[0]), 
				cheb::complex(bkr_y[0], bki_y[0])));
	}


cheb::Function globalFunction;
cheb::Function globalFunctionFirstDerivative;
cheb::Function globalFunctionSecondDerivative;

cheb::complex evalFunction(cheb::complex location){
    return clenshaw(location, globalFunction.funs[0].coeffs());
}
std::mutex m;
cheb::complex evalDerivative(cheb::complex location){
	auto [fx, Jx] = clenshawDeriv(location, globalFunction.funs[0].coeffs());
	auto c = clenshaw(location, globalFunctionFirstDerivative.funs[0].coeffs());

	//{
	//	auto fc = clenshaw(location, globalFunction.funs[0].coeffs());
	//	std::lock_guard lg(m);
	//	auto hx = 1e-8;
	//	auto x = clenshaw(location + cheb::complex(0., 0.), globalFunction.funs[0].coeffs());
	//	auto xr = clenshaw(location + cheb::complex(hx, 0.), globalFunction.funs[0].coeffs());
	//	auto xi = clenshaw(location + cheb::complex(0., hx), globalFunction.funs[0].coeffs());

	//	std::cout << location.real() << " : " << location.imag() << " -> " << fx.real() << " : " << fx.imag() << " - " << fc.real() << " - " << fc.imag() << std::endl;
	//	std::cout << "Clenshaw: " << c.real() << " : " << c.imag() << std::endl;
	//	std::cout << "Jacobian: " << Jx.first.real() << " : " << Jx.second.real() << " \\ " << Jx.first.imag() << " : " << Jx.second.imag() << std::endl;
	//	std::cout << "Finite:   " << (xr - x).real() / hx << " : " << (xr - x).imag() / hx << " \\ " << (xi - x).real() / hx << " : " << (xi - x).imag() / hx << std::endl;
	//	std::cout << std::endl;
	//	std::cout << c.real() / Jx.first.real() << " : " << c.imag() / Jx.second.real() << std::endl;
	//	std::cout << c.real() / Jx.first.imag() << " : " << c.imag() / Jx.second.imag() << std::endl;
 //
	//
	//
	//}
	return cheb::complex(Jx.first.real(), Jx.first.imag());
	return cheb::complex(Jx.second.imag(), -Jx.second.real());

    return clenshaw(location, globalFunctionFirstDerivative.funs[0].coeffs());
}
cheb::complex evalSecondDerivative(cheb::complex location){
    return clenshaw(location, globalFunctionSecondDerivative.funs[0].coeffs());
}