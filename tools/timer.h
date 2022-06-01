#pragma once
#include <algorithm>
#include <chrono>
// #if __has_include(<cuda_runtime.h>)
// #ifndef CUDA_SUPPORT
// #define CUDA_SUPPORT
// #endif
// #include <cuda_runtime.h>
// #endif
#include <map>
#include <mutex>
#include <random>
#include <set>
#include <string>
#include <boost/optional.hpp>
#include <tools/color.h>
#include <vector>
#include <numeric>
#include <array>

template<typename T = float, std::size_t bins = 32>
struct statistics {
    float avg, stddev;
    T min, max, sum, median;
    std::array<float, bins> histogram;

    statistics(std::vector<float> data) {
        if (data.size() < 2) return;
        std::sort(data.begin(), data.end());
        median = data[data.size() / 2];
        min = *data.begin();
        max = *data.rbegin();
        sum = std::accumulate(data.begin(), data.end(), 0.f);
        avg = (float)sum / (float)data.size();
        stddev = std::sqrt(std::accumulate(data.begin(), data.end(), 0.f, [&](float acc, T cur) {return acc + ((float)cur - avg) * ((float)cur - avg); }) / (float)data.size());
        for (auto& b : histogram)
            b = 0.f;
        for (auto e : data) {
            auto bin = std::round((float)(e - min) / (float)(max - min) * (float)(bins - 1));
            histogram[(std::size_t) bin]++;
        }
    }
};

class TimerManager;
class Timer {
public:
protected:
    std::vector<float> times;

    std::string identifier, sourceLocation = "";
    int32_t line = 0;
    Color color;

public:
    std::mutex vec_lock;
    bool graph;

    Timer(std::string descriptor = "unnamed", Color c = Color::black, std::string _location = "", int32_t _line = 0);
    virtual void start() = 0;
    virtual void stop() = 0;

    inline const std::vector<float>& getSamples() const { return times; }
    inline void setSamples(std::vector<float> arg) { times = arg; }
    inline std::string getDecriptor() { return identifier; }
    inline Color getColor() { return color; }
    inline std::string getSourceLocation() { return sourceLocation; }
    inline int32_t getLine() { return line; }

    boost::optional<statistics<float, 64>> getStats();
    friend class TImerManager;
};

class hostTimer : public Timer {
    std::chrono::high_resolution_clock::time_point last_tp;

public:
    hostTimer(std::string descriptor = "unnamed", Color c = Color::black, std::string _location = "", int32_t _line = 0);
    virtual void start() override;
    virtual void stop() override;
};
#ifdef CUDA_SUPPORT
class cudaTimer : public Timer {
    bool timer_stopped = false;
    cudaEvent_t start_event, stop_event;

public:
    cudaTimer(std::string descriptor = "unnamed", Color c = Color::black, std::string _location = "", int32_t _line = 0);
    virtual void start() override;
    virtual void stop() override;
};
extern std::map<std::string, cudaTimer*> timer_map;
#endif 

class TimerManager {
    static std::vector<Timer*> timers;

public:
    TimerManager() = delete;
    static hostTimer* createTimer(std::string name = "UNNAMED", Color c = Color::rosemadder,
        std::string _location = "", int32_t _line = 0);

    static Timer* createHybridTimer(std::string name = "UNNAMED", Color c = Color::rosemadder,
        std::string _location = "", int32_t _line = 0);
#ifdef CUDA_SUPPORT
    static cudaTimer* createGPUTimer(std::string name = "UNNAMED", Color c = Color::rosemadder,
        std::string _location = "", int32_t _line = 0);
#endif
    static inline std::vector<Timer*>& getTimers() { return timers; }
};


#define TIME_CODE(id, name, x)                                                                    \
  static hostTimer &t_## id = *TimerManager::createTimer(name, getHashedColor(name), __FILE__, __LINE__);                      \
  t_## id.start();                                                                                \
  x;                                                                                               \
  t_## id.stop();
#define GRAPH_CODE(name, col, x)                                                                   \
  static hostTimer &t_##name = *TimerManager::createTimer(#name, col,  __FILE__, __LINE__);                       \
  t_##name.start();                                                                                \
  x;                                                                                               \
  t_##name.stop();
#ifdef CUDA_SUPPORT
#define TIME_CODE_GPU(name, col, x)                                                                \
  static cudaTimer &t_##name = *TimerManager::createGPUTimer(#name, col,  __FILE__, __LINE__);                   \
  t_##name.start();                                                                                \
  x;                                                                                               \
  t_##name.stop();
#define GRAPH_CODE_GPU(name, col, x)                                                               \
  static cudaTimer &t_##name = *TimerManager::createGPUTimer(#name, col,  __FILE__, __LINE__);                    \
  t_##name.start();                                                                                \
  x;                                                                                               \
  t_##name.stop();
#endif