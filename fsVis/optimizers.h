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
#include <cheb/cheb.h>

using vec = Eigen::Vector2d;
using scalar = double;
using matrix = Eigen::Matrix2d;
using complex = std::complex<scalar>;
using clk = std::chrono::high_resolution_clock;

using v2 = Eigen::Vector2d;
using mat2 = Eigen::Matrix2d;

inline float *scalarFieldData = nullptr;
inline float *angularFieldData = nullptr;
inline float *normFieldData = nullptr;
inline int32_t dataWidth = 32;
inline int32_t dataHeight = 32;

extern cheb::Function globalFunction;

struct Jacobian {
  cheb::scalar dfdx = 0.0, dfdy = 0.0;
};

struct Hessian {
  cheb::scalar d2fdx2 = 0.0, d2fdxy = 0.0;
  cheb::scalar d2fdyx = 0.0, d2fdy2 = 0.0;
};

struct FunctionState {
  cheb::scalar f = 0.0;
  Jacobian J;
  Hessian H;
};
std::tuple<cheb::complex, cheb::complex, cheb::complex> evalPolynomial(cheb::complex location);
FunctionState evalSquarePolynomial(cheb::complex location);

enum struct optimizationMethod { newton, halley, gradientDescent, newtonOptimizer, newtonOptimizerHessian, Adam, gradientDescentAdaptive, adaGrad, BFGS };
enum optimizationState : int32_t { converged, diverged, unstable };

optimizationMethod getMethod(std::string stringMethod);
std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> BFGS(cheb::complex location);
std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> newtonsMethod(cheb::complex location);
std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> newtonsMethodOptimizer(cheb::complex location);
std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>>
newtonsMethodOptimizerHessian(cheb::complex location);
std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> gradientDescentAdaptive(cheb::complex location);
std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> adaGrad(cheb::complex location);
std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> gradientDescent(cheb::complex location);
std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> halleysMethod(cheb::complex location);
std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> adam(cheb::complex location);
std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> optimize(cheb::complex location,
                                                                                                                                          optimizationMethod method);