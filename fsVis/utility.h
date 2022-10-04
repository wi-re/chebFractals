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
#include <fsVis/utility.h>
#include <fsVis/optimizers.h>
#include <iomanip>
#include <cheb/cheb.h>

struct v4 {
  float x, y, z, w;
};

inline int32_t screenWidth = 2048 * std::sqrt(2);
inline int32_t screenHeight = 2048;
constexpr inline scalar epsilon((scalar)1e-7);

void initializeParameters();
void render();
void initRender();

uint32_t create1DTexture(v4 *colorMap, int32_t elements);
uint32_t createProgram(std::string vertexSource, std::string fragmentSource);
// converts a std::chrono duration to a millisecond value
scalar toMs(clk::duration dur);
std::pair<vec, vec> getDomain();

std::filesystem::path resolveFile(std::string fileName, std::vector<std::string> search_paths = {});

extern std::vector<std::vector<double>> globalCoefficients;
extern std::vector<std::tuple<cheb::complex, cheb::complex, cheb::complex>> trace;
