#include <algorithm>
#include <array>
#include <atomic>
#include <boost/range/combine.hpp>
#include <cfloat>
#include <chrono>
#include <fsVis/optimizers.h>
#include <fsVis/utility.h>
#include <iostream>
#include <iterator>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>

cheb::Function globalFunction;
std::vector<std::vector<double>> globalCoefficients;
scalar realOffset, complexOffset;

void initializeParameters() {
  static auto vectorModes = std::vector<detail::iAny>{std::string("real"),          std::string("imag"),         std::string("abs"),    std::string("real gradient"),
                                                      std::string("imag gradient"), std::string("abs gradient"), std::string("fractal")};
  static auto cMapPresets = std::vector<std::string>{[]() {
    std::vector<std::string> colorMaps;
    auto f = std::filesystem::path(ParameterManager::instance().get<std::string>("stylePath"));
    auto p = f.parent_path().string();
    if (*(p.end() - 1) == '/' || *(p.end() - 1) == '\\')
      p = p.substr(0, p.length() - 1);
    std::replace(p.begin(), p.end(), '\\', '/');
    for (auto &p : std::filesystem::directory_iterator(p))
      if (p.path().extension().string().find(".png") != std::string::npos)
        colorMaps.push_back(p.path().filename().replace_extension("").string());
    return colorMaps;
  }()};
  static auto colorMaps = std::vector<detail::iAny>{};
  for (auto c : cMapPresets)
    colorMaps.push_back(std::string(c));

  // {
  //         std::string("viridis"),
  //         std::string("hot"),
  //         std::string("RdBu") };

  ParameterManager::instance().newParameter("colorMap.renderMode", std::string("fractal"), {.constant = false, .presets = vectorModes});
  ParameterManager::instance().newParameter("colorMap.colorMap", std::string("twilight"), {.constant = false, .presets = colorMaps});

  auto genRange = [](auto min, auto max) { return Range{min, max}; };
  ParameterManager::instance().newParameter("index", 9, {.constant = false, .range = genRange(0, (int32_t)globalCoefficients.size() - 1)});

  ParameterManager::instance().newParameter("path.x", 0., {.constant = false, .range = genRange(-1., 1.)});
  ParameterManager::instance().newParameter("path.y", 0., {.constant = false, .range = genRange(-1., 1.)});
  ParameterManager::instance().newParameter("path.display", true, {.constant = false});

  // ParameterManager::instance().newParameter("domain.xmin", -1., { .constant = false ,.range = genRange(-10.,10.) });
  // ParameterManager::instance().newParameter("domain.xmax", 1., { .constant = false ,.range = genRange(-10.,10.) });
  ParameterManager::instance().newParameter("domain.xmin", -16. / 9., {.constant = false, .range = genRange(-10., 10.)});
  ParameterManager::instance().newParameter("domain.xmax", 16. / 9., {.constant = false, .range = genRange(-10., 10.)});
  ParameterManager::instance().newParameter("domain.ymin", -1., {.constant = false, .range = genRange(-10., 10.)});
  ParameterManager::instance().newParameter("domain.ymax", 1., {.constant = false, .range = genRange(-10., 10.)});

  // ParameterManager::instance().newParameter("domain.xmin", 1.5 * -16./9., { .constant = false ,.range = genRange(-10.,10.) });
  // ParameterManager::instance().newParameter("domain.xmax", 1.5 * 16./9., { .constant = false ,.range = genRange(-10.,10.) });
  // ParameterManager::instance().newParameter("domain.ymin", 1.5 * -1., { .constant = false ,.range = genRange(-10.,10.) });
  // ParameterManager::instance().newParameter("domain.ymax", 1.5 * 1., { .constant = false ,.range = genRange(-10.,10.) });

  ParameterManager::instance().newParameter("field.steps", 8192, {.constant = false, .range = genRange(0, 512)});

  ParameterManager::instance().newParameter("field.render", false, {.constant = false});
  ParameterManager::instance().newParameter("field.nx", 256, {.constant = false, .range = genRange(1, 512)});
  ParameterManager::instance().newParameter("field.ny", 256, {.constant = false, .range = genRange(1, 512)});
  ParameterManager::instance().newParameter("field.h", 1.0, {.constant = false, .range = genRange(0.01, 10.0)});
  ParameterManager::instance().newParameter("field.min", 0.0, {.constant = false, .range = genRange(-10.0, 10.0)});
  ParameterManager::instance().newParameter("field.max", 1.0, {.constant = false, .range = genRange(-10.0, 10.0)});
  ParameterManager::instance().newParameter("field.learningRate", -3.0, {.constant = false, .range = genRange(-8.0, 8.0)});

  ParameterManager::instance().newParameter("adam.alpha", .0, {.constant = false, .range = genRange(-8.0, 0.0)});
  ParameterManager::instance().newParameter("adam.beta1", 0.9, {.constant = false, .range = genRange(.0, 1.0)});
  ParameterManager::instance().newParameter("adam.beta2", 0.999, {.constant = false, .range = genRange(.0, 1.0)});
  ParameterManager::instance().newParameter("adam.eps", -8.0, {.constant = false, .range = genRange(-8.0, 0.0)});

  ParameterManager::instance().newParameter("BFGS.c1", 1e-4, {.constant = false, .range = genRange(.0, 1.0)});
  ParameterManager::instance().newParameter("BFGS.c2", 0.999, {.constant = false, .range = genRange(.0, 1.0)});
  ParameterManager::instance().newParameter("BFGS.alpha", 1., {.constant = false, .range = genRange(1e-4, 1.0)});
  ParameterManager::instance().newParameter("BFGS.eps", -8.0, {.constant = false, .range = genRange(-8.0, 0.0)});

  ParameterManager::instance().newParameter("colorMap.min", scalar(0.99), {.constant = false, .range = genRange(-10.0, 10.0)});
  ParameterManager::instance().newParameter("colorMap.max", scalar(1.03), {.constant = false, .range = genRange(-10.0, 10.0)});
  ParameterManager::instance().newParameter("colorMap.auto", true, {.constant = false});
  ParameterManager::instance().newParameter("field.clusterEpsilon", -2., {.constant = false, .range = genRange(-10., 1.)});
  ParameterManager::instance().newParameter("field.threshold", -2., {.constant = false, .range = genRange(-10., 1.)});
  ParameterManager::instance().newParameter("field.clustering", false, {.constant = false});
  ParameterManager::instance().newParameter("field.cycles", false, {.constant = false});
  static auto methods = std::vector<detail::iAny>{std::string("newton"), std::string("newton optimizer"), std::string("newton hessian"),
                                                  std::string("adam"),   std::string("halley"),           std::string("gradientDescent")};
  ParameterManager::instance().newParameter("field.method", std::string("newton"), {.constant = false, .presets = methods});
  ParameterManager::instance().newParameter("field.offset", scalar(0.00), { .constant = false, .range = genRange(-1.,1.)});
  ParameterManager::instance().newParameter("field.coffset", scalar(0.00), { .constant = false, .range = genRange(-1.,1.) });


  ParameterManager::instance().newParameter("recording.min_offset", scalar(-1.00), { .constant = false, .range = genRange(-1.,1.) });
  ParameterManager::instance().newParameter("recording.max_offset", scalar(1.00), { .constant = false, .range = genRange(-1.,1.) });
  ParameterManager::instance().newParameter("recording.step", scalar(0.01), { .constant = false, .range = genRange(-1.,1.) });
  ParameterManager::instance().newParameter("recording.active", bool(false), { .constant = false});
  ParameterManager::instance().newParameter("recording.done", bool(false), { .constant = false });

}

std::pair<vec, vec> getDomain() {
  static auto &xmin = ParameterManager::instance().get<scalar>("domain.xmin");
  static auto &xmax = ParameterManager::instance().get<scalar>("domain.xmax");
  static auto &ymin = ParameterManager::instance().get<scalar>("domain.ymin");
  static auto &ymax = ParameterManager::instance().get<scalar>("domain.ymax");
  return std::make_pair(vec(xmin, ymin), vec(xmax, ymax));
}

using clk = std::chrono::high_resolution_clock;
scalar toMs(clk::duration dur) { return static_cast<scalar>(std::chrono::duration_cast<std::chrono::microseconds>(dur).count()) / scalar(1000.0); }
std::vector<std::tuple<cheb::complex, cheb::complex, cheb::complex>> trace;
std::filesystem::path expand(std::filesystem::path in) {
  namespace fs = std::filesystem;
#ifndef _WIN32
  if (in.string().size() < 1)
    return in;

  const char *home = getenv("HOME");
  if (home == NULL) {
    std::cerr << "error: HOME variable not set." << std::endl;
    throw std::invalid_argument("error: HOME environment variable not set.");
  }

  std::string s = in.string();
  if (s[0] == '~') {
    s = std::string(home) + s.substr(1, s.size() - 1);
    return fs::path(s);
  } else {
    return in;
  }
#else
  if (in.string().size() < 1)
    return in;

  const char *home = getenv("USERPROFILE");
  if (home == NULL) {
    std::cerr << "error: USERPROFILE variable not set." << std::endl;
    throw std::invalid_argument("error: USERPROFILE environment variable not set.");
  }

  std::string s = in.string();
  if (s[0] == '~') {
    s = std::string(home) + s.substr(1, s.size() - 1);
    return fs::path(s);
  } else {
    return in;
  }
#endif
}
std::filesystem::path resolveFile(std::string fileName, std::vector<std::string> search_paths) {
  namespace fs = std::filesystem;
  auto &pm = ParameterManager::instance();
  fs::path working_dir = pm.get<std::string>("internal.working_directory");
  fs::path binary_dir = pm.get<std::string>("internal.binary_directory");
  fs::path source_dir = pm.get<std::string>("internal.source_directory");
  fs::path build_dir = pm.get<std::string>("internal.build_directory");
  fs::path expanded = expand(fs::path(fileName));

  fs::path base_path = "";
  if (fs::exists(expand(fs::path(fileName))))
    return expand(fs::path(fileName));
  for (const auto &path : search_paths) {
    auto p = expand(fs::path(path));
    if (fs::exists(p / fileName))
      return p.string() + std::string("/") + fileName;
  }

  if (fs::exists(fileName))
    return fs::path(fileName);
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
