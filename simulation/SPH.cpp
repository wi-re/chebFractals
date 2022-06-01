
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
            bk1 = ak[0] + x * bk1 - bk2;
            bk2 = temp;
		}
		cheb::complex res = 0.0;
		res = ak[0] + 0.5 * x * bk1 - bk2;

		return res;
	}


cheb::Function globalFunction;

cheb::complex evalFunction(cheb::complex location){
    static bool once = true;
    if(once){
        once = false;
    }

    cheb::svec coeffs = globalFunction.funs[0].coeffs();
    return clenshaw(location, coeffs);
}