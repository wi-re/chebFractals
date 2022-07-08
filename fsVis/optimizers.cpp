#include <fsVis/utility.h>
#include <fsVis/optimizers.h>
#include <algorithm>
#include <array>
#include <boost/range/combine.hpp>
#include <iostream>
#include <numeric>
#include <chrono>
#include <sstream>
#include <atomic>
#include <random>
#include <iterator>
#include <cfloat>
#include <mutex>

cheb::complex clenshaw(cheb::complex x, const cheb::svec &ak) {
  if (ak.size() == 0)
    return 0.0;
  if (ak.size() == 1)
    return ak[0];
  // if (any(ak, [](auto s) {return std::isnan(s); }))
  //     return std::nan("") * svec(1., xx.size());
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

struct complexState {
  FunctionState r, i;
};

complexState clenshawDeriv(cheb::complex val, const cheb::svec &ak) {
  if (ak.size() == 0)
    return complexState{};
  if (ak.size() == 1)
    return complexState{.r = FunctionState{.f = ak[0]}};

  std::vector<complexState> b(ak.size() + 2);
  auto x = val.real();
  auto y = val.imag();

  for (int32_t k = (int32_t)ak.size() - 1; k >= 0; k -= 1) {
    auto s = k == 0 ? 1. : 2.;

    b[k].r.f = ak[k] - b[k + 2].r.f + s * (x * b[k + 1].r.f - y * b[k + 1].i.f);

    b[k].r.J.dfdx = -b[k + 2].r.J.dfdx + s * (x * b[k + 1].r.J.dfdx + b[k + 1].r.f - y * b[k + 1].i.J.dfdx);

    b[k].r.H.d2fdx2 = -b[k + 2].r.H.d2fdx2 + s * (x * b[k + 1].r.H.d2fdx2 + b[k + 1].r.J.dfdx + b[k + 1].r.J.dfdx - y * b[k + 1].i.H.d2fdx2);
    b[k].r.H.d2fdxy = -b[k + 2].r.H.d2fdxy + s * (x * b[k + 1].r.H.d2fdxy + b[k + 1].r.J.dfdy - y * b[k + 1].i.H.d2fdxy - b[k + 1].i.J.dfdx);

    b[k].r.J.dfdy = -b[k + 2].r.J.dfdy + s * (x * b[k + 1].r.J.dfdy - y * b[k + 1].i.J.dfdy - b[k + 1].i.f);

    b[k].r.H.d2fdyx = -b[k + 2].r.H.d2fdyx + s * (x * b[k + 1].r.H.d2fdyx + b[k + 1].r.J.dfdy - y * b[k + 1].i.H.d2fdyx - b[k + 1].i.J.dfdx);
    b[k].r.H.d2fdy2 = -b[k + 2].r.H.d2fdy2 + s * (x * b[k + 1].r.H.d2fdy2 - y * b[k + 1].i.H.d2fdy2 - b[k + 1].i.J.dfdy - b[k + 1].i.J.dfdy);

    b[k].i.f = -b[k + 2].i.f + s * (x * b[k + 1].i.f + y * b[k + 1].r.f);

    b[k].i.J.dfdx = -b[k + 2].i.J.dfdx + s * (x * b[k + 1].i.J.dfdx + b[k + 1].i.f + y * b[k + 1].r.J.dfdx);

    b[k].i.H.d2fdx2 = -b[k + 2].i.H.d2fdx2 + s * (x * b[k + 1].i.H.d2fdx2 + b[k + 1].i.J.dfdx + b[k + 1].i.J.dfdx + y * b[k + 1].r.H.d2fdx2);
    b[k].i.H.d2fdxy = -b[k + 2].i.H.d2fdxy + s * (x * b[k + 1].i.H.d2fdxy + b[k + 1].i.J.dfdy + y * b[k + 1].r.H.d2fdxy + b[k + 1].r.J.dfdx);

    b[k].i.J.dfdy = -b[k + 2].i.J.dfdy + s * (x * b[k + 1].i.J.dfdy + y * b[k + 1].r.J.dfdy + b[k + 1].r.f);

    b[k].i.H.d2fdyx = -b[k + 2].i.H.d2fdyx + s * (x * b[k + 1].i.H.d2fdyx + b[k + 1].i.J.dfdy + y * b[k + 1].r.H.d2fdyx + b[k + 1].r.J.dfdx);
    b[k].i.H.d2fdy2 = -b[k + 2].i.H.d2fdy2 + s * (x * b[k + 1].i.H.d2fdy2 + y * b[k + 1].r.H.d2fdy2 + b[k + 1].r.J.dfdy + b[k + 1].r.J.dfdy);
  }

  return b[0];
}

std::tuple<cheb::complex, cheb::complex, cheb::complex> evalPolynomial(cheb::complex location) {
  auto [fr, fi] = clenshawDeriv(location, globalFunction.funs[0].coeffs());

  cheb::complex f(fr.f, fi.f);
  Jacobian J(fr.J.dfdx, fi.J.dfdx);
  Hessian H(fr.H.d2fdx2, fr.H.d2fdxy, fi.H.d2fdx2, fi.H.d2fdxy);
  return std::make_tuple(f, cheb::complex(J.dfdx, J.dfdy), cheb::complex(H.d2fdx2, H.d2fdyx));
}

FunctionState evalSquarePolynomial(cheb::complex location) {
  auto [fr, fi] = clenshawDeriv(location, globalFunction.funs[0].coeffs());

  auto f = fr.f * fr.f + fi.f * fi.f;

  auto dfdx = 2. * fr.f * fr.J.dfdx + 2. * fi.f * fi.J.dfdx;
  auto dfdy = 2. * fr.f * fr.J.dfdy + 2. * fi.f * fi.J.dfdy;

  auto d2fdx2 = 2. * (fr.J.dfdx * fr.J.dfdx + fr.f * fr.H.d2fdx2 + fi.J.dfdx * fi.J.dfdx + fi.f * fi.H.d2fdx2);
  auto d2fdxy = 2. * (fr.J.dfdy * fr.J.dfdx + fr.f * fr.H.d2fdxy + fi.J.dfdy * fi.J.dfdx + fi.f * fi.H.d2fdxy);
  auto d2fdyx = 2. * (fr.J.dfdy * fr.J.dfdx + fr.f * fr.H.d2fdyx + fi.J.dfdy * fi.J.dfdx + fi.f * fi.H.d2fdyx);
  auto d2fdy2 = 2. * (fr.J.dfdy * fr.J.dfdy + fr.f * fr.H.d2fdy2 + fi.J.dfdy * fi.J.dfdy + fi.f * fi.H.d2fdy2);

  Jacobian J{dfdx, dfdy};
  Hessian H{d2fdx2, d2fdxy, d2fdyx, d2fdy2};

  return FunctionState{f, J, H};
}

optimizationMethod getMethod(std::string stringMethod) {
  if (stringMethod == "newton")
    return optimizationMethod::newton;
  if (stringMethod == "halley")
    return optimizationMethod::halley;
  if (stringMethod == "gradientDescent")
    return optimizationMethod::gradientDescent;
  if (stringMethod == "gradientDescentAdaptive")
    return optimizationMethod::gradientDescentAdaptive;
  if (stringMethod == "newton optimizer")
    return optimizationMethod::newtonOptimizer;
  if (stringMethod == "newton hessian")
    return optimizationMethod::newtonOptimizerHessian;
  if (stringMethod == "adam")
    return optimizationMethod::Adam;
  if (stringMethod == "adaGrad")
    return optimizationMethod::adaGrad;
  if (stringMethod == "BFGS")
    return optimizationMethod::BFGS;
  return optimizationMethod::newton;
}
auto evalSQFunction(Eigen::Vector2d location) {
  auto [f, J, H] = evalSquarePolynomial(cheb::complex(location.x(), location.y()));
  return f;
}

// BFGS

auto finiteGrad(Eigen::Vector2d p) {
  auto [f, J, H] = evalSquarePolynomial(cheb::complex(p.x(), p.y()));
  return v2(J.dfdx, J.dfdy);

  auto h = std::cbrt(DBL_EPSILON);

  auto dx = (evalSQFunction(v2{p.x() + h, p.y()})) - evalSQFunction(v2{p.x() - h, p.y()}) / (2. * h);
  auto dy = (evalSQFunction(v2{p.x(), p.y() + h})) - evalSQFunction(v2{p.x(), p.y() - h}) / (2. * h);

  return v2(dx, dy);
}

auto lineSearchInv(Eigen::Vector2d x, Eigen::Vector2d p) {
  auto f = [](Eigen::Vector2d x) { return evalSQFunction(x); };
  auto nabla = [](Eigen::Vector2d x) { return finiteGrad(x); };

  auto alpha = ParameterManager::instance().get<scalar>("BFGS.alpha");
  static auto &c1 = ParameterManager::instance().get<scalar>("BFGS.c1");
  static auto &c2 = ParameterManager::instance().get<scalar>("BFGS.c2");
  static auto &epsp = ParameterManager::instance().get<scalar>("BFGS.eps");
  auto eps = std::pow(10., epsp);
  // auto alpha = -1.;
  // auto c1 = 1e-4;
  // auto c2 = 0.9;
  int32_t i = 0;

  bool armijo = f(x + alpha * p) <= f(x) + c1 * alpha * p.dot(nabla(x));
  bool curvature = -p.dot(nabla(x + alpha * p)) <= -c2 * p.dot(nabla(x));

  while ((!armijo || !curvature) && std::abs(alpha) >= eps) {
    alpha *= 0.5;
    auto x_new = x + alpha * p;
    auto nabla_new = finiteGrad(x_new);
    armijo = f(x + alpha * p) <= f(x) + c1 * alpha * p.dot(nabla(x));
    curvature = -p.dot(nabla(x + alpha * p)) <= -c2 * p.dot(nabla(x));
    // printf(
    //     "\t[%03d]: [%g %g] + %g * [%g %g] = [%g %g] -> [%g %g] -> %g >= %g || %g <= %g\n",
    //     ++i, x.x(), x.y(), alpha, p.x(), p.y(), x_new.x(), x_new.y(), nabla_new.x(), nabla_new.y(),
    //     f(x + alpha * p) , f(x) + c1 * alpha * p.dot(nabla(x)),
    //     -p.dot(nabla(x + alpha * p)) , - c2 * p.dot(nabla(x)));
  }
  return alpha;
}
auto lineSearch(Eigen::Vector2d x, Eigen::Vector2d p) {
  auto f = [](Eigen::Vector2d x) { return evalSQFunction(x); };
  auto nabla = [](Eigen::Vector2d x) { return finiteGrad(x); };
  auto alpha = ParameterManager::instance().get<scalar>("BFGS.alpha");
  static auto &c1 = ParameterManager::instance().get<scalar>("BFGS.c1");
  static auto &c2 = ParameterManager::instance().get<scalar>("BFGS.c2");
  static auto &epsp = ParameterManager::instance().get<scalar>("BFGS.eps");
  auto eps = std::pow(10., epsp);

  int32_t i = 0;

  bool armijo = f(x + alpha * p) <= f(x) + c1 * alpha * p.dot(nabla(x));
  bool curvature = -p.dot(nabla(x + alpha * p)) <= -c2 * p.dot(nabla(x));

  while ((!armijo || !curvature) && alpha >= eps) {
    alpha *= 0.5;
    auto x_new = x + alpha * p;
    auto nabla_new = finiteGrad(x_new);
    armijo = f(x + alpha * p) <= f(x) + c1 * alpha * p.dot(nabla(x));
    curvature = -p.dot(nabla(x + alpha * p)) <= -c2 * p.dot(nabla(x));
    // printf(
    //     "\t[%03d]: [%g %g] + %g * [%g %g] = [%g %g] -> [%g %g] -> %g >= %g || %g <= %g\n",
    //     ++i, x.x(), x.y(), alpha, p.x(), p.y(), x_new.x(), x_new.y(), nabla_new.x(), nabla_new.y(),
    //     f(x + alpha * p) , f(x) + c1 * alpha * p.dot(nabla(x)),
    //     -p.dot(nabla(x + alpha * p)) , - c2 * p.dot(nabla(x)));
  }
  if (alpha < eps)
    return lineSearchInv(x, p);
  return alpha;
}

std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> BFGS(cheb::complex location) {
  std::vector<cheb::complex> positions, values, steps;
  static auto &stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
  static auto &learningRatef = ParameterManager::instance().get<scalar>("field.learningRate");

  auto x0 = v2{location.real(), location.imag()};

  auto [f, J, H] = evalSquarePolynomial(cheb::complex(location.real(), location.imag()));

  mat2 Bk;
  // Bk << H.d2fdx2, H.d2fdxy, H.d2fdyx, H.d2fdy2;
  Bk << 1., 0., 0., 1.;
  Bk = Bk.inverse();
  v2 x = x0;
  int32_t it = 2;

  positions.push_back(cheb::complex(x.x(), x.y()));
  values.push_back(evalSQFunction(x));
  auto step = cheb::complex(0., 0.);
  steps.push_back(step);

  for (int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i) {
    it += 1;
    v2 p = -Bk * finiteGrad(x);
    // printf("[%03d]: %g + %gi => Starting Linesearch\n",it, p.x(), p.y());

    auto a = lineSearch(x, p);
    auto s = a * p;
    // printf("[%03d]: a = %g : [%g + %gi]\n",it, a, s.x(), s.y());
    v2 x_new = x + s;

    auto y = finiteGrad(x_new) - finiteGrad(x);

    scalar sTy = s.dot(y);
    Eigen::RowVector2d yTB = y.transpose() * Bk;
    scalar yTBy = yTB * y;
    mat2 ssT = s * s.transpose();
    mat2 leftTerm = (sTy + yTBy) / (sTy * sTy) * ssT;

    mat2 BysT = (Bk * y) * s.transpose();
    mat2 syTB = (s * y.transpose()) * Bk;
    mat2 rightTerm = (BysT + syTB) / (sTy);

    // mat2 Bkp = Bk + leftTerm - rightTerm;

    mat2 Id = mat2::Identity();
    scalar rhok = 1. / y.dot(s);
    mat2 left = Id - rhok * s * y.transpose();
    mat2 right = Id - rhok * y * s.transpose();
    mat2 righter = rhok * s * s.transpose();
    mat2 Bkp = left * Bk * right + righter;

    // printf("[%03d]: f(%g + %gi) = %g : [%g  + %gi] -> [%g %g] [%g %g]\n",it, x_new.x(), x_new.y(),evalSQFunction(x_new), finiteGrad(x_new).x(), finiteGrad(x_new).y(), Bkp(0,0),
    // Bkp(0,1), Bkp(1,0), Bkp(1,1));

    // mat2 Bkp = Bk + (s.transpose() * y + y.transpose() * Bk * y) *(s * s.transpose()) / (s.transpose() * y).squaredNorm() - (Bk * y * s.transpose() + s * y.transpose() * Bk) /
    // (s.transpose() * y)[0];

    positions.push_back(cheb::complex(x.x(), x.y()));
    values.push_back(evalSQFunction(x));
    auto step = cheb::complex(s.x(), s.y());
    steps.push_back(step);

    if (std::abs(step) < 1e-14)
      return std::make_tuple(optimizationState::converged, location, positions, values, steps);
    if (location != location || step != step)
      return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
    location += step;
    x = x_new;

    Bk = Bkp;
  }
  return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}

std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> newtonsMethod(cheb::complex location) {
  std::vector<cheb::complex> positions, values, steps;
  static auto &stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
  static auto &learningRatef = ParameterManager::instance().get<scalar>("field.learningRate");
  auto learningRate = std::pow(10., learningRatef);

  for (int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i) {
    auto [fx, dx, H] = evalPolynomial(location);
    auto step = -fx / dx;
    positions.push_back(location);
    values.push_back(fx);
    steps.push_back(step);

    if (std::abs(step) < 1e-12)
      return std::make_tuple(optimizationState::converged, location, positions, values, steps);
    if (location != location || step != step)
      return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
    location += step;
  }
  return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}
std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> newtonsMethodOptimizer(cheb::complex location) {
  std::vector<cheb::complex> positions, values, steps;
  static auto &stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
  static auto &learningRatef = ParameterManager::instance().get<scalar>("field.learningRate");
  auto learningRate = std::pow(10., learningRatef);

  for (int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i) {
    auto [f, J, H] = evalSquarePolynomial(location);

    std::complex a = H.d2fdx2;
    std::complex b = H.d2fdxy;
    std::complex c = H.d2fdyx;
    std::complex d = H.d2fdy2;
    std::complex T = a + d;
    std::complex D = a * d - b * c;
    std::complex l1 = T / 2. + std::sqrt(T * T / 4. - D);
    std::complex l2 = T / 2. - std::sqrt(T * T / 4. - D);
    auto det = H.d2fdx2 * H.d2fdy2 - H.d2fdxy * H.d2fdyx;

    Hessian H_1{H.d2fdy2 / det, -H.d2fdxy / det, -H.d2fdyx / det, H.d2fdx2 / det};

    auto prodx = H_1.d2fdx2 * J.dfdx + H_1.d2fdxy * J.dfdy;
    auto prody = H_1.d2fdyx * J.dfdx + H_1.d2fdy2 * J.dfdy;

    auto dx = cheb::complex(prodx, prody);
    auto step = -dx * 1.;

    positions.push_back(location);
    values.push_back(f);
    steps.push_back(step);
    if (std::abs(step) < 1e-12)
      return std::make_tuple(optimizationState::converged, location, positions, values, steps);
    if (location != location || step != step)
      return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
    location += step;
  }
  return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}
std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>>
newtonsMethodOptimizerHessian(cheb::complex location) {
  std::vector<cheb::complex> positions, values, steps;
  static auto &stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
  static auto &learningRatef = ParameterManager::instance().get<scalar>("field.learningRate");
  auto learningRate = std::pow(10., learningRatef);

  for (int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i) {
    auto [f, J, H] = evalSquarePolynomial(location);

    std::complex a = H.d2fdx2;
    std::complex b = H.d2fdxy;
    std::complex c = H.d2fdyx;
    std::complex d = H.d2fdy2;
    std::complex T = a + d;
    std::complex D = a * d - b * c;
    std::complex l1 = T / 2. + std::sqrt(T * T / 4. - D);
    std::complex l2 = T / 2. - std::sqrt(T * T / 4. - D);
    auto l = std::min(l1.real(), l2.real());
    if (l < 0.) {
      H.d2fdx2 -= l * 1.5;
      H.d2fdy2 -= l * 1.5;
    }

    auto det = H.d2fdx2 * H.d2fdy2 - H.d2fdxy * H.d2fdyx;

    Hessian H_1{H.d2fdy2 / det, -H.d2fdxy / det, -H.d2fdyx / det, H.d2fdx2 / det};

    auto prodx = H_1.d2fdx2 * J.dfdx + H_1.d2fdxy * J.dfdy;
    auto prody = H_1.d2fdyx * J.dfdx + H_1.d2fdy2 * J.dfdy;

    auto dx = cheb::complex(prodx, prody);
    auto step = -dx * 1.;

    positions.push_back(location);
    values.push_back(f);
    steps.push_back(step);
    if (std::abs(step) < 1e-12)
      return std::make_tuple(optimizationState::converged, location, positions, values, steps);
    if (location != location || step != step)
      return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
    location += step;
  }
  return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}
std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> gradientDescentAdaptive(cheb::complex location) {
  std::vector<cheb::complex> positions, values, steps;
  static auto &stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
  static auto &learningRatef = ParameterManager::instance().get<scalar>("field.learningRate");
  auto learningRate = std::pow(10., learningRatef);

  auto [f, J, H] = evalSquarePolynomial(location);
  auto dxPrior = cheb::complex(J.dfdx, J.dfdy);
  auto locationPrior = location;

  auto step = -dxPrior * learningRate;
  positions.push_back(location);
  values.push_back(f);
  steps.push_back(step);

  location += step;

  for (int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i) {
    auto [f, J, H] = evalSquarePolynomial(location);

    auto dx = cheb::complex(J.dfdx, J.dfdy);
    // dx = cheb::complex(prodx, prody);

    auto diff = location - locationPrior;
    auto diffGrad = dx - dxPrior;

    auto gamma = (diff.real() * diffGrad.real() + diff.imag() * diffGrad.imag()) / (diffGrad.real() * diffGrad.real() + diffGrad.imag() * diffGrad.imag() + 1e-12);

    auto step = -gamma * dx;

    positions.push_back(location);
    values.push_back(f);
    steps.push_back(step);
    if (std::abs(step) < 1e-16)
      return std::make_tuple(optimizationState::converged, location, positions, values, steps);
    if (location != location || step != step)
      return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
    locationPrior = location;
    location += step;
    dxPrior = dx;
  }
  return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}

std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> adaGrad(cheb::complex location) {
  std::vector<cheb::complex> positions, values, steps;
  static auto &stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
  static auto &learningRatef = ParameterManager::instance().get<scalar>("field.learningRate");
  auto learningRate = std::pow(10., learningRatef);

  auto sumx = 0.;
  auto sumy = 0.;

  for (int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i) {
    auto [f, J, H] = evalSquarePolynomial(location);

    auto dx = cheb::complex(J.dfdx, J.dfdy);
    // dx = cheb::complex(prodx, prody);

    sumx += J.dfdx * J.dfdx;
    sumy += J.dfdy * J.dfdy;

    auto alphax = learningRate / (1e-8 + std::sqrt(sumx));
    auto alphay = learningRate / (1e-8 + std::sqrt(sumy));
    auto step = -cheb::complex(J.dfdx * alphax, J.dfdy * alphay);

    // auto step = -dx * learningRate;

    positions.push_back(location);
    values.push_back(f);
    steps.push_back(step);
    if (std::abs(step) < 1e-12)
      return std::make_tuple(optimizationState::converged, location, positions, values, steps);
    if (location != location || step != step)
      return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
    location += step;
  }
  return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}
std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> gradientDescent(cheb::complex location) {
  std::vector<cheb::complex> positions, values, steps;
  static auto &stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
  static auto &learningRatef = ParameterManager::instance().get<scalar>("field.learningRate");
  auto learningRate = std::pow(10., learningRatef);

  for (int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i) {
    auto [f, J, H] = evalSquarePolynomial(location);

    auto dx = cheb::complex(J.dfdx, J.dfdy);
    // dx = cheb::complex(prodx, prody);

    auto step = -dx * learningRate;

    positions.push_back(location);
    values.push_back(f);
    steps.push_back(step);
    if (std::abs(step) < 1e-12)
      return std::make_tuple(optimizationState::converged, location, positions, values, steps);
    if (location != location || step != step)
      return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
    location += step;
  }
  return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}
std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> halleysMethod(cheb::complex location) {
  std::vector<cheb::complex> positions, values, steps;
  static auto &stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
  static auto &learningRatef = ParameterManager::instance().get<scalar>("field.learningRate");
  auto learningRate = std::pow(10., learningRatef);

  for (int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i) {
    auto [f, J, H] = evalPolynomial(location);
    auto step = -(2. * f * J) / (2. * J * J - f * H);

    positions.push_back(location);
    values.push_back(f);
    steps.push_back(step);

    if (std::abs(step) < 1e-12)
      return std::make_tuple(optimizationState::converged, location, positions, values, steps);
    if (location != location || step != step)
      return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
    location += step;
  }
  return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}

std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> adam(cheb::complex location) {
  std::vector<cheb::complex> positions, values, steps;
  static auto &stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
  static auto &alphaP = ParameterManager::instance().get<scalar>("adam.alpha");
  static auto &beta1 = ParameterManager::instance().get<scalar>("adam.beta1");
  static auto &beta2 = ParameterManager::instance().get<scalar>("adam.beta2");
  static auto &epsP = ParameterManager::instance().get<scalar>("adam.eps");

  auto alpha = std::pow(10., alphaP);
  auto eps = std::pow(10., epsP);

  auto t = 0.0;
  auto m_x = 0.;
  auto v_x = 0.;

  auto m_y = 0.;
  auto v_y = 0.;

  for (int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i) {
    auto [f, J, H] = evalSquarePolynomial(location);

    t = t + 1.;
    // get gradients
    auto g_t = J;

    m_x = beta1 * m_x + (1. - beta1) * g_t.dfdx;
    v_x = beta2 * v_x + (1. - beta2) * g_t.dfdx * g_t.dfdx;
    auto m_x_corrected = m_x / (1. - std::pow(beta1, t));
    auto v_x_corrected = v_x / (1. - std::pow(beta2, t));
    auto step_x = -alpha * m_x_corrected / (std::sqrt(v_x_corrected) + eps);

    m_y = beta1 * m_y + (1. - beta1) * g_t.dfdy;
    v_y = beta2 * v_y + (1. - beta2) * g_t.dfdy * g_t.dfdy;
    auto m_y_corrected = m_y / (1. - std::pow(beta1, t));
    auto v_y_corrected = v_y / (1. - std::pow(beta2, t));
    auto step_y = -alpha * m_y_corrected / (std::sqrt(v_y_corrected) + eps);

    auto step = cheb::complex(step_x, step_y);

    positions.push_back(location);
    values.push_back(f);
    steps.push_back(step);
    if (std::abs(step) < 1e-12)
      return std::make_tuple(optimizationState::converged, location, positions, values, steps);
    if (location != location || step != step)
      return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
    location += step;
  }
  return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}
std::tuple<optimizationState, cheb::complex, std::vector<cheb::complex>, std::vector<cheb::complex>, std::vector<cheb::complex>> optimize(cheb::complex location,
                                                                                                                                          optimizationMethod method) {
  switch (method) {
  case optimizationMethod::newton:
    return newtonsMethod(location);
  case optimizationMethod::newtonOptimizer:
    return newtonsMethodOptimizer(location);
  case optimizationMethod::newtonOptimizerHessian:
    return newtonsMethodOptimizerHessian(location);
  case optimizationMethod::gradientDescent:
    return gradientDescent(location);
  case optimizationMethod::gradientDescentAdaptive:
    return gradientDescentAdaptive(location);
  case optimizationMethod::halley:
    return halleysMethod(location);
  case optimizationMethod::adaGrad:
    return adaGrad(location);
  case optimizationMethod::Adam:
    return adam(location);
  case optimizationMethod::BFGS:
    return BFGS(location);
  default:
    return newtonsMethod(location);
  }
}