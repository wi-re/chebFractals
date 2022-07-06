#pragma once
//#define _USE_MATH_DEFINES
#include <array>
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <tools/ParameterManager.h>
#include <iostream>
#include <math.h>
#include <sstream>
#include <vector>
#include <fstream>
#include <numbers>

// comment this line out to use double precision for everything in the simulation
//#define USE_SINGLE
#ifndef USE_SINGLE
using vec = Eigen::Vector2d;
using scalar = double;
using matrix = Eigen::Matrix2d;
using complex = std::complex<scalar>;
#else
using vec = Eigen::Vector2f;
using scalar = float;
using matrix = Eigen::Matrix2f;
using complex = std::complex<scalar>;
#endif
// using define for easier usage of time
using clk = std::chrono::high_resolution_clock;
constexpr inline scalar domainScale = 20.0;

inline std::ofstream summaryFile;
inline bool summaryFileOpen = false;

void initializeParameters(int32_t scene = 0);
// execute a single SPH timestep
void render();
void initRender();
// converts a std::chrono duration to a millisecond value
scalar toMs(clk::duration dur);

// inline int32_t screenWidth = 2048;
// inline int32_t screenHeight = 2048;

inline int32_t screenWidth = 3840;
inline int32_t screenHeight = 2160;


inline float* scalarFieldData = nullptr;
inline float* angularFieldData = nullptr;
inline float* normFieldData = nullptr;
inline int32_t dataWidth = 32;
inline int32_t dataHeight = 32;


// used to create a gap from the domain to the edge of the window

std::pair<vec,vec> getDomain();
#include <iomanip>


// numerical Parameters
constexpr inline scalar epsilon((scalar)1e-7);


struct Triangle {
    vec v0 = vec(0, 0), v1 = vec(0, 0), v2 = vec(0, 0);
};


extern  std::vector<std::vector<double>> globalCoefficients;

std::filesystem::path resolveFile(std::string fileName, std::vector<std::string> search_paths = {});
#include <cheb/cheb.h>
extern std::vector<std::tuple<cheb::complex, cheb::complex, cheb::complex>> trace;

extern cheb::Function globalFunction;
extern cheb::Function globalFunctionFirstDerivative;
extern cheb::Function globalFunctionSecondDerivative;


	struct Jacobian{
		cheb::scalar dfdx = 0.0, dfdy = 0.0; 
	};

	struct Hessian{
		cheb::scalar d2fdx2 = 0.0, d2fdxy = 0.0;
		cheb::scalar d2fdyx = 0.0, d2fdy2 = 0.0;
	};

	struct FunctionState{
		cheb::scalar f = 0.0;
		Jacobian J;
		Hessian H;
	};

cheb::complex evalFunction(cheb::complex location);
cheb::complex evalDerivative(cheb::complex location);
cheb::complex evalSecondDerivative(cheb::complex location);

std::tuple<cheb::complex, cheb::complex, cheb::complex>  evalPolynomial(cheb::complex location);
FunctionState evalSquarePolynomial(cheb::complex location);