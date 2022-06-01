#include <tools/exceptionHandling.h>
#include "gui/glui.h"
#include <sstream>
#include <thread>
//#include <windows.h>
#include <iostream> 

// BOOL CtrlHandler(DWORD fdwCtrlType) {
//     std::clog << "Caught signal " << fdwCtrlType << std::endl;
//     switch (fdwCtrlType)
//     {
//     case CTRL_CLOSE_EVENT:
//         GUI::instance().quit();
//         GUI::instance().render_lock.unlock();
//         std::this_thread::sleep_for(std::chrono::milliseconds(2000));
//         return(TRUE);

//     default:
//         return FALSE;
//     }
// }
#include <test.h>
#include <simulation/2DMath.h>
#include <cheb/cheb.h>
#include <iomanip>

#include <tools/timer.h>
int main(int argc, char* argv[]) 
{
    globalFunction = cheb::Function([](cheb::scalar x){ return x * x * x;},cheb::Domain{-1.,1.});

    auto& gui = GUI::instance();
    gui.render_lock.lock();
    gui.initParameters(argc, argv);
    gui.initGVDB();
    gui.initSimulation();
    gui.initGL();
    gui.renderLoop();


    for (auto tptr : TimerManager::getTimers()) {
        auto& t = *tptr;
        auto stats = t.getStats().value();
        scalar sum = 0.;
        for (auto s : t.getSamples()) {
            sum += s;
        }
        std::cout << std::setprecision(4);
        std::cout << t.getDecriptor() << ": " << stats.median << "ms / " << stats.avg << "ms @ " << stats.stddev << " : " << stats.min << "ms / " << stats.max << "ms, total: " << sum/1000. << "s\n";

    }

    return 0;
}