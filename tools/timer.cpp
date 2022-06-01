#include <algorithm>
#include <tools/timer.h>
#ifdef CUDA_SUPPORT
std::map<std::string, cudaTimer*> timer_map;
cudaTimer* TimerManager::createGPUTimer(std::string name, Color c, std::string s, int32_t l) {
    if (c == Color::rosemadder)
        c = getRandomColor();
    cudaTimer* new_timer = new cudaTimer(name, c, s, l);
    timers.push_back(new_timer);
    return new_timer;
}

Timer* TimerManager::createHybridTimer(std::string name, Color c, std::string s, int32_t l) {
    if (get<parameters::internal::target>() == launch_config::device)
        return createGPUTimer(name, c, s, l);
    else
        return createTimer(name, c, s, l);
}


cudaTimer::cudaTimer(std::string descr, Color c, std::string s, int32_t l) : Timer(descr, c, s, l) {
    cudaEventCreate(&start_event);
    cudaEventCreate(&stop_event);
}

void cudaTimer::start() {
    if (timer_stopped) {
        cudaEventSynchronize(stop_event);
        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start_event, stop_event);

        std::lock_guard<std::mutex> guard(vec_lock);
        times.push_back((milliseconds));
        timer_stopped = false;
    }
    cudaEventRecord(start_event);
}
void cudaTimer::stop() {
    cudaEventRecord(stop_event);

    timer_stopped = true;
}

#endif

boost::optional<statistics<float, 64>> Timer::getStats() {
    if (times.size() == 0) return boost::none;
    return statistics<float, 64>(times);
}

Timer::Timer(std::string descr, Color c, std::string loc, int32_t l)
    : identifier(descr), color(c), sourceLocation(loc), line(l) {}

std::vector<Timer*> TimerManager::timers;

hostTimer::hostTimer(std::string descr, Color c, std::string s, int32_t l) : Timer(descr, c, s, l) {
    last_tp = std::chrono::high_resolution_clock::now();
}

hostTimer* TimerManager::createTimer(std::string name, Color c, std::string s, int32_t l) {
    if (c == Color::rosemadder)
        c = getRandomColor();
    hostTimer* new_timer = new hostTimer(name, c, s, l);
    timers.push_back(new_timer);
    return new_timer;
}
void hostTimer::start() { last_tp = std::chrono::high_resolution_clock::now(); }
void hostTimer::stop() {
    auto tp = std::chrono::high_resolution_clock::now();
    auto dur_us = std::chrono::duration_cast<std::chrono::microseconds>(tp - last_tp).count();
    std::lock_guard<std::mutex> guard(vec_lock);
    times.push_back((static_cast<float>(dur_us) / 1000.f));
}
