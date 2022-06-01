#include <simulation/SPH.h>
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
#include <random>
#include <iterator>
#include <cfloat>
#include <mutex>

void initializeParameters(int32_t scene) {
    static auto vectorModes = std::vector<detail::iAny>{
            std::string("real"),
            std::string("imag"),
            std::string("real gradient"),
            std::string("imag gradient"),
            std::string("newton") };
    static auto cMapPresets = std::vector<std::string>{ []() {std::vector <std::string> colorMaps; auto f = std::filesystem::path(ParameterManager::instance().get<std::string>("stylePath")); auto p = f.parent_path().string(); if (*(p.end() - 1) == '/' || *(p.end() - 1) == '\\')p = p.substr(0, p.length() - 1); std::replace(p.begin(), p.end(), '\\', '/'); for (auto& p : std::filesystem::directory_iterator(p))if (p.path().extension().string().find(".png") != std::string::npos)colorMaps.push_back(p.path().filename().replace_extension("").string()); return colorMaps; }() };
    static auto colorMaps = std::vector<detail::iAny>{};
    for(auto c : cMapPresets)
    colorMaps.push_back(std::string(c));
    
    // {
    //         std::string("viridis"),
    //         std::string("hot"),
    //         std::string("RdBu") };

    ParameterManager::instance().newParameter("colorMap.renderMode", std::string("magnitude"), { .constant = false, .presets =vectorModes
        });
    ParameterManager::instance().newParameter("colorMap.colorMap", std::string("viridis"), { .constant = false, .presets = colorMaps
        });

    ParameterManager::instance().newParameter("domain.xmin", -1., { .constant = false ,.range = Range{-10.,10.}});
    ParameterManager::instance().newParameter("domain.xmax", 1., { .constant = false ,.range = Range{-10.,10.}});
    ParameterManager::instance().newParameter("domain.ymin", -1., { .constant = false ,.range = Range{-10.,10.}});
    ParameterManager::instance().newParameter("domain.ymax", 1., { .constant = false ,.range = Range{-10.,10.}});

    ParameterManager::instance().newParameter("field.render", false, { .constant = false });
    ParameterManager::instance().newParameter("field.nx", 32, { .constant = false , .range = Range{1,512} });
    ParameterManager::instance().newParameter("field.ny", 32, { .constant = false, .range = Range{1,512}  });
    ParameterManager::instance().newParameter("field.h", 1.0, { .constant = false, .range = Range{0.01,10.0}  });
    ParameterManager::instance().newParameter("field.min", 0.0, { .constant = false , .range = Range{-10.0,10.0} });
    ParameterManager::instance().newParameter("field.max", 1.0, { .constant = false, .range = Range{-10.0,10.0} });
    
    ParameterManager::instance().newParameter("colorMap.min", scalar(0.99), { .constant = false , .range = Range{-10.0,10.0} });
    ParameterManager::instance().newParameter("colorMap.max", scalar(1.03), { .constant = false , .range = Range{-10.0,10.0} });
    ParameterManager::instance().newParameter("colorMap.auto", true, { .constant = false });

}


std::pair<vec,vec> getDomain(){
    static auto& xmin = ParameterManager::instance().get<scalar>("domain.xmin");
    static auto& xmax = ParameterManager::instance().get<scalar>("domain.xmax");
    static auto& ymin = ParameterManager::instance().get<scalar>("domain.ymin");
    static auto& ymax = ParameterManager::instance().get<scalar>("domain.ymax");
    return std::make_pair(vec(xmin,ymin), vec(xmax,ymax));
}