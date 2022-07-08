#include "gui/glui.h"
#include <sstream>
#include <thread>
#include <iostream> 
#include <cheb/cheb.h>
#include <iomanip>
#include <vector>
#include <fstream>
#include <filesystem>


int main(int argc, char* argv[]) 
{
    // f(x) = cos(13x) * sin(25x)
    globalFunction = cheb::Function([](cheb::scalar x){ return std::cos(x * 13.) * sin(x * 25.);},cheb::Domain{-1.,1.});
    globalCoefficients.push_back(globalFunction.funs[0].coeffs());
    // f(x) = tan(-x)
    globalFunction = cheb::Function([](cheb::scalar x) { return std::tan(x); }, cheb::Domain{-1.,1.});
    globalCoefficients.push_back(globalFunction.funs[0].coeffs());
    // f(x) = abs(x)
    // cheb::eps = std::numeric_limits<scalar>::epsilon() * 1e3;
    // globalFunction = cheb::Function([](cheb::scalar x) { return std::abs(x); }, cheb::Domain{-1.,1.});
    // globalCoefficients.push_back(globalFunction.funs[0].coeffs());
    // cheb::eps = std::numeric_limits<scalar>::epsilon();
    // f(x) = exp(-x)sin(x)
    globalFunction = cheb::Function([](cheb::scalar x) { return std::exp(-x) * std::sin(x); }, cheb::Domain{-1.,1.});
    globalCoefficients.push_back(globalFunction.funs[0].coeffs());
    // f(x) = x^3 - 2x + 1
    globalFunction = cheb::Function([](cheb::scalar x) { return x * x * x - 2. * x + 2.; }, cheb::Domain{ -1.,1. });
    globalCoefficients.push_back(globalFunction.funs[0].coeffs());
    // f(x) = x^3 - 1
    globalFunction = cheb::Function([](cheb::scalar x) { return x * x * x - 1.; }, cheb::Domain{ -1.,1. });
    globalCoefficients.push_back(globalFunction.funs[0].coeffs());
    // f(x) = x^5 - 1
    globalFunction = cheb::Function([](cheb::scalar x){ return x * x * x * x * x - .5; }, cheb::Domain{-1., 1.});
    globalCoefficients.push_back(globalFunction.funs[0].coeffs());
    // f(x) = tanh(x)
    globalFunction = cheb::Function([](cheb::scalar x){ return std::tanh(x); }, cheb::Domain{-1., 1.});
    globalCoefficients.push_back(globalFunction.funs[0].coeffs());
    // f(x) = ln(1 + exp(x))
    globalFunction = cheb::Function([](cheb::scalar x){ return std::log(1 + std::exp(x)); }, cheb::Domain{-1., 1.});
    globalCoefficients.push_back(globalFunction.funs[0].coeffs());
    // f(x) = 0.5 x (1 + erf(x/sqrt(2)))
    globalFunction = cheb::Function([](cheb::scalar x){ return 0.5 * x * (1 + std::erf(x / std::sqrt(2))); }, cheb::Domain{-1., 1.});
    globalCoefficients.push_back(globalFunction.funs[0].coeffs());
    // f(x) = erf(x)
    globalFunction = cheb::Function([](cheb::scalar x){ return std::erf(x); }, cheb::Domain{-1., 1.});
    globalCoefficients.push_back(globalFunction.funs[0].coeffs());
    // Teaser image
    globalCoefficients.push_back({0x1.cac0684e08fa9p-7, 0x1.175db08d85ee6p-29, 0x1.7f1a95ebf3356p-35, -0x1.4a3a6350cbe08p-37, 0x1.4938bc58e7dcap-46, 0x1.23653faff6428p-45, -0x1.33226145f4p-51, -0x1.bccc0a6fp-54 });
    // Test function 717
    globalCoefficients.push_back({0x1.10fc7c5569b9ep+3, -0x1.51437cca43795p+0, -0x1.f5bab4a18965ep-2, -0x1.562b1df168a7cp-3, -0x1.c20092cd47efcp-5, -0x1.2d2e61ba9f892p-6, -0x1.99b9c624554a8p-8, -0x1.1a65df8552b52p-9, -0x1.8973970604da4p-11, -0x1.147bb89c7a32dp-12, -0x1.8761f07a244f4p-14, -0x1.16b145ee47474p-15, -0x1.8ef041975995ep-17, -0x1.1ecb5ac5f861p-18, -0x1.9dee1a4615aap-20, -0x1.2bb72ff2c0586p-21, -0x1.b351a054c6ab3p-23, -0x1.3cf890831d49ap-24, -0x1.ceaf60800749p-26, -0x1.52690ac65ec6p-27, -0x1.effb86808cbp-29, -0x1.6c190ff6861p-30, -0x1.0bb57658bap-31, -0x1.8a4090078p-33, -0x1.22b229bep-34, -0x1.ad38812p-36, -0x1.3d3b0ep-37, -0x1.d57352p-39, -0x1.5bc06p-40, -0x1.01bf8p-41, -0x1.7dc5p-43, -0x1.1ab8p-44, -0x1.ap-46});
    // Test function 7173
    globalCoefficients.push_back({0x1.501dd4b3eac76p+3, -0x1.1b2382387516ep+0, 0x1.063f1e052e0bep-1, -0x1.354fbf4835015p-7, -0x1.a79e986a0a751p-9, -0x1.cc1c4b94b001ap-11, -0x1.d4ed4411b1f4dp-13, -0x1.d2b50f32433adp-15, -0x1.cbf6735b6be15p-17, -0x1.c382ee845d5c4p-19, -0x1.baa5d1e76f96p-21, -0x1.b1faebe41b78p-23, -0x1.a9ca3aa428p-25, -0x1.a23219a71p-27, -0x1.9b3bce08p-29, -0x1.94e5eaap-31, -0x1.8f298p-33, -0x1.89fdc6p-35, -0x1.8559cp-37, -0x1.814076p-39, -0x1.7da262p-41, -0x1.7b96b98p-43, -0x1.740529p-45});
    // Test function 71734
    globalCoefficients.push_back({0x1.1ab97ffe88f6dp+3, -0x1.1eff319b19f4ep+0, -0x1.27c1920bb214fp-1, -0x1.fd4654d5e066dp-3, -0x1.8de1b703fbeb5p-4, -0x1.3ab9cd1da7a16p-5, -0x1.f8ca1677a832ap-7, -0x1.99a612f24be6cp-8, -0x1.4fba5c7a24f55p-9, -0x1.156608e715a8p-10, -0x1.cd8a3357dff0dp-12, -0x1.822c90756f89ep-13, -0x1.44b4ef3e7a6f4p-14, -0x1.122efb11a112p-15, -0x1.d0c1b3ebe51fep-17, -0x1.8b2cba84d4983p-18, -0x1.50fa14025671bp-19, -0x1.20161b9c215a7p-20, -0x1.edb5d18cdb788p-22, -0x1.a7ed2452ef68p-23, -0x1.6cb08244488cp-24, -0x1.3a4426dbf8f8p-25, -0x1.0f3cdf753fcp-26, -0x1.d4df6b1188p-28, -0x1.95ca9c46cp-29, -0x1.5fa01128p-30, -0x1.3108cc08p-31, -0x1.08e50f1p-32, -0x1.cc86584p-34, -0x1.90ae72p-35, -0x1.5ceb2p-36, -0x1.301504p-37, -0x1.094p-38, -0x1.cee28p-40, -0x1.943ccp-41, -0x1.61688p-42, -0x1.35158p-43, -0x1.0db38p-44, -0x1.dbcp-46});
    // Test function 717346
    globalCoefficients.push_back({0x1.e1a9bab7850ccp+1, -0x1.632ff31acd9e7p+0, 0x1.a51b400959a36p+1, -0x1.108a2013afb52p-1, -0x1.6af899ca4697dp-2, -0x1.d1727786ab0a5p-3, -0x1.2857c8a431afep-3, -0x1.79c7ee5732c44p-4, -0x1.e34a8fc627522p-5, -0x1.366572c3f8139p-5, -0x1.904f93669fc86p-6, -0x1.031a9991e2f04p-6, -0x1.508fa45f4dfe4p-7, -0x1.b686caab210a2p-8, -0x1.1e7d5dd597f1dp-8, -0x1.774598f3f0d04p-9, -0x1.ecb06160ceb12p-10, -0x1.4417805e2f641p-10, -0x1.ab2e1ef5d1324p-11, -0x1.1a03a715408fcp-11, -0x1.74f215e50625cp-12, -0x1.ede958921102fp-13, -0x1.477e5b840c92ap-13, -0x1.b2d5509f16dfap-14, -0x1.2102104efbc92p-14, -0x1.809514a3645e3p-15, -0x1.0022a1c482cbp-15, -0x1.557e2e786c68fp-16, -0x1.c7b0e23cc858ap-17, -0x1.304942ad1ef8ep-17, -0x1.96aeccc1cae06p-18, -0x1.0ff6d75c6acc5p-18, -0x1.6bfe0d4df5f11p-19, -0x1.e7793cffc22e5p-20, -0x1.469eca02a19fap-20, -0x1.b5f0a33b44ecep-21, -0x1.25c29328a20fp-21, -0x1.8a4c8ec990488p-22, -0x1.08c10a23a84ap-22, -0x1.63b52b15a726p-23, -0x1.de1f1148ec8ep-24, -0x1.4177f2ced4cap-24, -0x1.b075c9882ac8p-25, -0x1.230006a11e74p-25, -0x1.87c56fc98b2p-26, -0x1.07d055c974cp-26, -0x1.636b7dd828p-27, -0x1.defe37b11cp-28, -0x1.42de159cfp-28, -0x1.b364b451ap-29, -0x1.25a794cfp-29, -0x1.8c3a27b08p-30, -0x1.0b631055p-30, -0x1.68fa3a76p-31, -0x1.e7731a3cp-32, -0x1.49328118p-32, -0x1.bcbfcb4p-33, -0x1.2c7fd1ep-33, -0x1.9628e6p-34, -0x1.128bbbcp-34, -0x1.733bdp-35, -0x1.f614dcp-36, -0x1.53966p-36, -0x1.cb7544p-37, -0x1.36dep-37, -0x1.a4be9p-38, -0x1.1ccacp-38, -0x1.818ccp-39, -0x1.050ecp-39, -0x1.6186p-40, -0x1.dee1ap-41, -0x1.44612p-41, -0x1.b769cp-42, -0x1.29b78p-42, -0x1.93864p-43, -0x1.12112p-43, -0x1.739aep-44, -0x1.f8578p-45, -0x1.555a8p-45, -0x1.d1782p-46, -0x1.3cb2ep-46, -0x1.b02d4p-47, -0x1.2fe7dp-47});
    // Test function 110109
    globalCoefficients.push_back({0x1.2d68d860a1dfap+3, -0x1.d2fcf66941e4ap-1, -0x1.26094d5cca5ecp-1, -0x1.03a4799511c02p-2, -0x1.9957f994885f9p-4, -0x1.469660f91d278p-5, -0x1.08209d8a36004p-6, -0x1.b04a08e216b38p-8, -0x1.653ce3a2ef52p-9, -0x1.299f7d9c1b3bcp-10, -0x1.f34a698d89781p-12, -0x1.a536869b92ep-13, -0x1.65173d65b7adap-14, -0x1.3003729395a3ap-15, -0x1.03c7fd9ecf651p-16, -0x1.bd68b898b7187p-18, -0x1.7eef6f914bef7p-19, -0x1.4a1105d690058p-20, -0x1.1d264716c4352p-21, -0x1.edb5915d27abp-23, -0x1.ac34dbc8b524p-24, -0x1.7407320cf22p-25, -0x1.43b96e7b966p-26, -0x1.1a187f905p-27, -0x1.ec49ed5d3p-29, -0x1.ae12ba688p-30, -0x1.782516ccp-31, -0x1.49532eep-32, -0x1.209d4c2p-33, -0x1.fa5656p-35, -0x1.bc892ep-36, -0x1.869a68p-37, -0x1.578p-38, -0x1.2e387p-39, -0x1.0a246p-40, -0x1.d4828p-42, -0x1.9d76p-43, -0x1.6b498p-44, -0x1.3c858p-45});
    // Test function 1487127
    globalCoefficients.push_back({0x1.abc42004e3953p+0, 0x1.bbb01fbcaefa2p+1, 0x1.082a9fd68def4p+2, 0x1.738ef1c2e6747p-13, -0x1.a146d43740c54p-16, -0x1.2bef0c399955ep-18, -0x1.28b4ace759101p-21, -0x1.07f708eadab4p-24, -0x1.c133d6aee3199p-28, -0x1.75db6b631131p-31, -0x1.336eb9aa5f8p-34, -0x1.f625322538p-38, -0x1.98210188p-41, -0x1.48e4a18p-44});
    // Test function 340364
    globalCoefficients.push_back({0x1.f01aa39e426a3p+3, 0x1.3ff214b721929p+0, 0x1.ac64c6cd3146ep-8, -0x1.35c6645efb989p-7, -0x1.92dc57b26c46p-10, -0x1.0148c36d813d8p-12, -0x1.48d7f200a2592p-15, -0x1.a65e3f47e7236p-18, -0x1.10c693bc82e38p-20, -0x1.623ac47d5548p-23, -0x1.ce4079fbbd5p-26, -0x1.2eead8bf1p-28, -0x1.8e88d28d4p-31, -0x1.070f8c84p-33, -0x1.5c53d5p-36, -0x1.ce968p-39, -0x1.34p-41, -0x1.a114p-44});

    // std::cout << "Coefficient count: " << globalFunction.funs[0].coeffs().size() << std::endl;

    auto& gui = GUI::instance();
    gui.render_lock.lock();
    gui.initParameters(argc, argv);
    gui.initGVDB();
    gui.initSimulation();
    gui.initGL();
    gui.renderLoop();


    return 0;
}