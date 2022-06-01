#pragma once
//#define _USE_MATH_DEFINES
#include "SPH.h"
#include <math.h>
// this set of power functions is used to avoid any potential inaccuracies of calling a normal pow function with integer
// powers and always uses a simple iterative multiplication of the base.
inline auto power(std::complex<scalar> v, int32_t p) {
    if (p == 0)
        return std::complex<scalar>(scalar(0.0), scalar(0.0));
    std::complex<scalar> s = v;
    for (int32_t i = 1; i < p; ++i)
        s *= v;
    return s;
}
inline auto power(scalar v, int32_t p) {
    if (p == 0)
        return scalar(0.0);
    scalar s = v;
    for (int32_t i = 1; i < p; ++i)
        s *= v;
    return s;
}
// for some functions it is easier to describe complex values as 2.0 + 1.0_i instead of requiring std::complex(2,1)
inline constexpr std::complex<scalar> operator"" _i(long double d) {
    return std::complex<scalar>(scalar(0.0), static_cast<scalar>(d));
}
inline constexpr std::complex<scalar> operator"" _i(unsigned long long d) {
    return std::complex<scalar>(scalar(0.0), static_cast<scalar>(d));
}
// checks if a given scalar is zero, i.e. abs(scalar) < epsilon
inline auto isZero(scalar v) { return std::abs(v) < epsilon; }
// returns the sign of the given scalar. 1 if v > 0, -1 if v < 1.0, 0 else
inline scalar sgn(scalar v) { return v < scalar(-epsilon) ? scalar(-1.0) : (v > scalar(epsilon) ? scalar(1.0) : scalar(0.0)); }
// computes the SVD of a given 2x2 matrix using a closed form solution
inline std::tuple<matrix, matrix, matrix> svd2x2(matrix M) {
    scalar a = M(0, 0);
    scalar a2 = a * a;
    scalar b = M(0, 1);
    scalar b2 = b * b;
    scalar c = M(1, 0);
    scalar c2 = c * c;
    scalar d = M(1, 1);
    scalar d2 = d * d;

    scalar theta = scalar(0.5) * std::atan2(scalar(2.0) * a * c + scalar(2.0) * b * d, a2 + b2 - c2 - d2);
    scalar ct = std::cos(theta);
    scalar st = std::sin(theta);
    matrix U;
    U << ct, -st, st, ct;

    scalar S1 = a2 + b2 + c2 + d2;
    scalar S2 = std::sqrt(power(a2 + b2 - c2 - d2, 2) + scalar(4.0) * power(a * c + b * d, 2));
    scalar sigma1 = std::sqrt((S1 + S2) * scalar(0.5));
    scalar sigma2 = std::sqrt((S1 - S2) * scalar(0.5));
    matrix S;
    S << sigma1, scalar(0.0), scalar(0.0), sigma2;

    scalar phi = scalar(0.5) * std::atan2(scalar(2.0) * a * b + scalar(2.0) * c * d, a2 - b2 + c2 - d2);
    auto cp = std::cos(phi);
    auto sp = std::sin(phi);
    scalar s11 = (a * ct + c * st) * cp + (b * ct + d * st) * sp;
    scalar s22 = (a * st - c * ct) * sp + (-b * st + d * ct) * cp;
    matrix V;
    V << sgn(s11) * cp, -sgn(s22) * sp, sgn(s11)* sp, sgn(s22)* cp;

    return std::make_tuple(U, S, V);
}