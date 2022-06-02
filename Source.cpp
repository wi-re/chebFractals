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
    std::vector<scalar> coefficients{ 0x1.0000000000000p+0, -0x1.bfb5dbce01c41p-5, 0x1.08a2242d1a359p-1, 0x1.d7d54c363507dp-4, -0x1.5ab0400be2843p-2, -0x1.d1eae3a019217p-4, 0x1.53467431323d9p-3, 0x1.5eca3b803d81fp-4, -0x1.58e81cada5e99p-5, -0x1.e6236ce4c4606p-5, -0x1.c12efc8eac5e8p-6, 0x1.6276d0bf5035ep-5, 0x1.bec289393e0e2p-5, -0x1.1a219e6b86f2cp-5, -0x1.be6f9dca8d40bp-5, 0x1.c86edb95cfa5ep-6, 0x1.4faaec94f4b4cp-5, -0x1.531681d75560ep-6, -0x1.8b071512834abp-6, 0x1.9eec819b561b8p-7, 0x1.4685577420891p-7, -0x1.407ccf235ffb4p-8, -0x1.f826949ba693ap-11, -0x1.4d20028064b9ep-10, -0x1.8f467b5f5dfb1p-9, 0x1.4848b898bdbdcp-8, 0x1.cf8dedbe3c789p-9, -0x1.9b8c05d32b86cp-8, -0x1.2a668ad17ef8fp-9, 0x1.722b800f2f559p-8, 0x1.6163325b83c0dp-11, -0x1.03deb36dc62fdp-8, 0x1.ed102d5493f9cp-12, 0x1.0c0a899298079p-9, -0x1.e0e11a665aa61p-11, -0x1.e1fcde06475d1p-12, 0x1.a17d424cfcff5p-11, -0x1.1420db519d807p-11, -0x1.9620e087d4793p-12, 0x1.debd026d02dd2p-11, -0x1.77e2d9b1cfb88p-15, -0x1.c2a4c616293c8p-11, 0x1.5fa4e33273f1cp-12, 0x1.2f25b30b64a58p-11, -0x1.c728c92780f56p-12, -0x1.09bebd3bcdd60p-12, 0x1.8969dcd9e5390p-12, 0x1.680273ede1100p-20, -0x1.e695886d22a44p-13, 0x1.19d3db31d1c22p-13, 0x1.451d862a51092p-14, -0x1.5830b483f469ep-13, 0x1.47feae7218d71p-15, 0x1.080519441cabcp-13, -0x1.9dbd1095e7c3dp-14, -0x1.040581ab6cdc2p-14, 0x1.bdd7093daa787p-14, 0x1.1edc0b9d64c65p-18, -0x1.51cc6693f858ep-14, 0x1.05395b4abaeebp-15, 0x1.6338e2adf00abp-15, -0x1.613a8e4c9e76cp-15, -0x1.0a23005ae65a4p-17, 0x1.2713dca1a13fcp-15, -0x1.dc2544c9d7998p-17, -0x1.4d0ded3ce256ep-16, 0x1.7b9b2a8f7e71bp-16, 0x1.2590172e562fap-18, -0x1.59f326c541e35p-16, 0x1.a9b4323a1c73fp-18, 0x1.ba5470220865dp-17, -0x1.6d1c5c6a91d44p-17, -0x1.463a52fd827b3p-18, 0x1.5a32cadecb44bp-17, -0x1.80b6e24fbbfefp-20, -0x1.cd39e56e96a2fp-18, 0x1.366773fa23429p-18, 0x1.71cab65c2fec4p-19, -0x1.5205f7fca7e6fp-18, 0x1.14c353fe0d9d4p-21, 0x1.f16b4fcc3e700p-19, -0x1.33f1dd603abf4p-19, -0x1.dbb6cf7610516p-20, 0x1.64084d89122d9p-19, 0x1.6dd3f47632fabp-24, -0x1.150322e20639cp-19, 0x1.f5c7d0c7726b0p-21, 0x1.26036598cdf96p-20, -0x1.4ec6fff1a3e10p-20, -0x1.a61dad27fe1dap-23, 0x1.190d666bced26p-20, -0x1.a1469b963d869p-22, -0x1.45a66f0c91cf0p-21, 0x1.48ad8898fca18p-21, 0x1.5796bcdd57b4cp-23, -0x1.2a75cb7178f83p-21, 0x1.489cc531d2b5cp-23, 0x1.7af0d6daebba4p-22, -0x1.38b650ef06663p-22, -0x1.0e96f5eaffd60p-23, 0x1.30eb66ea77323p-22, -0x1.86dd2765c282fp-25, -0x1.9eb96c1b94f0dp-23, 0x1.1b87f08109b25p-23, 0x1.5ab30fd3869f0p-24, -0x1.316dec25e342cp-23, 0x1.59ac252061764p-27, 0x1.c0f26c684ecb2p-24, -0x1.037615fb53c96p-24, -0x1.aff1985341829p-25, 0x1.3385568a55678p-24, 0x1.975029396f7d1p-29, -0x1.e3f215e378686p-25, 0x1.bada768c2a27dp-26, 0x1.03418d1154959p-25, -0x1.2acce3e83b422p-25, -0x1.812e5674328bep-28, 0x1.f9930aa1acfeep-26, -0x1.6a0b052bdd143p-27, -0x1.289a48b8d3084p-26, 0x1.207e243827d7fp-26, 0x1.4b580bdc3afa5p-28, -0x1.075a94e0fe19ap-26, 0x1.113de10ebf786p-28, 0x1.5006ec8f4406dp-27, -0x1.10b92ed181238p-27, -0x1.e4cdf13420658p-29, 0x1.0d6b308d1b188p-27, -0x1.50acc4addcea3p-30, -0x1.723cbcdc933dcp-28, 0x1.f48281f359180p-29, 0x1.3bd4524d638a6p-29, -0x1.0f3992e52404dp-28, 0x1.e66b0c53c1961p-33, 0x1.910480e16b5f9p-29, -0x1.c152ed8ced550p-30, -0x1.863b4865a9378p-30, 0x1.0e2d141dd36aap-29, 0x1.cb1835183e3ccp-34, -0x1.aca6f786691b9p-30 };


    // globalFunction = cheb::Function([](cheb::scalar x){ return std::cos(x * 13.) * sin(x * 25.);},cheb::Domain{-1.,1.});
    // globalFunction = cheb::Function(std::vector<cheb::IntervalFunction>{cheb::IntervalFunction(coefficients)});
    // globalFunction = cheb::Function([](cheb::scalar x) { return std::tan(x); }, cheb::Domain{-1.,1.});
    globalFunction = cheb::Function([](cheb::scalar x) { return x * x * x - 1.; }, cheb::Domain{ -1.,1. });
    // globalFunction = cheb::Function([](cheb::scalar x) { return x * x * x - 2. * x + 2.; }, cheb::Domain{ -1.,1. });
    globalFunctionFirstDerivative = globalFunction.derivative();
    globalFunctionSecondDerivative = globalFunctionFirstDerivative.derivative();

    std::cout << evalFunction(cheb::complex(1., 0.)) << " :8 " << globalFunction(1.) << std::endl;

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