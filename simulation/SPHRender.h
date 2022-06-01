#pragma once
#include <glad/glad.h>
#include <GLFW/glfw3.h>
//#define _USE_MATH_DEFINES
#include <array>
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <tools/ParameterManager.h>
#include <iostream>
#include <math.h>
#include <sstream>
#include <vector>
#include <glm/glm.hpp>

struct v4 { float x, y, z, w; };
GLuint create1DTexture(v4* colorMap, int32_t elements);
GLuint createProgram(std::string vertexSource, std::string fragmentSource);

void renderMarching();
void renderRays();


inline auto square(scalar x) { return x * x; }
inline auto cubic(scalar x) { return x * x * x; }
inline auto renderKernel(scalar q) {
	return std::max(0.0, cubic(1.0 - q * q));
}
inline auto positionKernel(scalar q) {
	return std::max(0.0, 1.0 - cubic(q));
}
inline std::size_t cellsMX, cellsMY;
#include <simulation/2DMath.h>