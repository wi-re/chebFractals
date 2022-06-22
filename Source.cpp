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


struct progressBar {
    std::chrono::high_resolution_clock::time_point start, lastTime;
    int32_t frames, lastFrame = 0;

    progressBar(int32_t targetFrames) :frames(targetFrames) {
        start = std::chrono::high_resolution_clock::now();
    };

    void update(int32_t current, int32_t updateRate) {
        std::ios cout_state(nullptr);
        cout_state.copyfmt(std::cout);
        if (current >= lastFrame + updateRate) {
            auto progress = ((double)current) / ((double)frames);
            auto now = std::chrono::high_resolution_clock::now();
            lastTime = now;
            int barWidth = 70;
            std::cout << "Rendering " << std::setw(4) << current;
            if (frames != -1)
                std::cout << "/" << std::setw(4) << frames;
            std::cout << " [";
            int pos = barWidth * progress;
            for (int i = 0; i < barWidth; ++i) {
                if (i < pos)
                    std::cout << "=";
                else if (i == pos)
                    std::cout << ">";
                else
                    std::cout << " ";
            }
            std::cout << "] " << std::setw(3) << int(progress * 100.0) << " ";
            auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
            if (dur.count() < 100 || progress < 1e-3f) {
                std::cout << " ---/---s  ";
            }
            else {
                auto totalTime =
                    ((float)std::chrono::duration_cast<std::chrono::microseconds>(now - start).count()) / 1000.f / 1000.f;
                std::cout << std::fixed << std::setprecision(0) << " " << std::setw(3) << totalTime << "/" << std::setw(3)
                    << (totalTime / progress) << "s  ";
            }


        }
        std::cout << "\r";
        std::cout.flush();
        std::cout.copyfmt(cout_state);
    }
    void update(int32_t x, int32_t y, int32_t updateRate) {
        std::ios cout_state(nullptr);
        cout_state.copyfmt(std::cout);
        int32_t current = x + y * 1920;
        if (current >= lastFrame + updateRate) {
            auto progress = ((double)current) / ((double)frames);
            auto now = std::chrono::high_resolution_clock::now();
            lastTime = now;
            int barWidth = 70;
            std::cout << "Rendering " << std::setw(5) << x << std::setw(5) << y;
            if (frames != -1)
                std::cout << "/" << std::setw(5) << 1920 << std::setw(5) << 1080;
            std::cout << " [";
            int pos = barWidth * progress;
            for (int i = 0; i < barWidth; ++i) {
                if (i < pos)
                    std::cout << "=";
                else if (i == pos)
                    std::cout << ">";
                else
                    std::cout << " ";
            }
            std::cout << "] " << std::setw(3) << int(progress * 100.0) << " ";
            auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
            if (dur.count() < 100 || progress < 1e-3f) {
                std::cout << " ---/---s  ";
            }
            else {
                auto totalTime =
                    ((float)std::chrono::duration_cast<std::chrono::microseconds>(now - start).count()) / 1000.f / 1000.f;
                std::cout << std::fixed << std::setprecision(0) << " " << std::setw(3) << totalTime << "/" << std::setw(3)
                    << (totalTime / progress) << "s  ";
            }


        }
        std::cout << "\r";
        std::cout.flush();
        std::cout.copyfmt(cout_state);
    }
    void end() {
        std::cout << std::endl;
    }
};
#include <vector>
#include <fstream>
#include <filesystem>

std::vector<std::vector<double>> chebyshevData(std::string file) {
    
constexpr std::size_t width = 1920, height = 1080;

struct float4 { float x, y, z, w; };
struct double3 { double x, y, z; };
struct float3 { float x, y, z; };
struct int3 { int32_t x, y, z; };

struct Matrix4x4d {
    double data[4 * 4] = {
        0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0. };
};

struct isoFunctionState {
    // iso function information
    double(*isoFunctionA)(double, double, double);  // function pointer to first iso function (A)
    double(*isoFunctionB)(double, double, double);  // function pointer to second iso function (B)
    Matrix4x4d matrixA =                            // transformation matrix applied to isofunction A
    { 1.,0.,0.,0.,
     0.,1.,0.,0.,
     0.,0.,1.,0.,
     0.,0.,0.,1. };
    Matrix4x4d matrixB =                            // transformation matrix applied to isofunction B
    { 1.,0.,0.,0.,
     0.,1.,0.,0.,
     0.,0.,1.,0.,
     0.,0.,0.,1. };
    // blending and isosurface parameters
    // The next parameter describes the blend mode to use when evaluating the scalar field with:
    // blendmode =  0 : A
    // blendmode =  1 : B
    // blendmode =  2 : A + B
    // blendmode =  3 : A - B
    // blendmode =  4 : A * B + blendBias
    // blendmode =  5 : 0.5 ( A + B - sqrt(A^2 + B^2))
    // blendmode =  6 : 0.5 ( A + B + sqrt(A^2 + B^2))
    // blendmode =  7 : A / B
    // blendmode =  8 : blendParameter * A + (1 - blendParameter) * B
    // blendmode =  9 : max(A,B)
    // blendmode = 10 : min(A,B)
    // blendmode = 11 : - blendA0 / (1 + A^2 / blendA1 + B^2 / blendA2)
    int32_t blendMode = 0;
    double blendParameter = 0.5;                    // linear blend weight for blendmode 8
    double blendBias = 0.01;                        // bias parameter for multiplying isosurfaces for blendmode 4
    double blendA0 = 0.019;                         // parameter for blendmode 11
    double blendA1 = 0.589;                         // parameter for blendmode 11
    double blendA2 = 0.943;                         // parameter for blendmode 11
    double isoLevel = 0.0;                          // offset value applied to any evaluated isofunction, set to negative for blendmode 11

    // gradient settings
    int32_t diffScheme = 0;                         // flag to pick gradient scheme: 0 is central finite, 1 is forward finite, 3 is chebyshev
    int32_t order = 2;                              // order for finite difference schemes (2,4,6,8 for central; 1,2,3,4,5,6 for forward)
    bool autoScale = true;                          // automatically scale step width for x as h = sqrt(x_+ - x)
    double h = DBL_EPSILON;                         // manual step width for finite differences
    // sturm sequence settings
    int32_t splitDegree = 50;                       // degree at which proxy functions should be split for computational efficiency
    bool useQSeries = false;                        // utilize quotient based sturm sequences instead of remainder based
    double thresholdA = 12.0;                       // numerical precision threshold during sturm sequence construction
    double thresholdB = 500.0;                      // numerical precision threshold during sturm sequence construction
    // domain settings
    float3 minDomain = float3{ -1.5f, -1.5f, -1.5f }; // minimum of AABB surrounding the isofunction
    float3 maxDomain = float3{ 1.5f,  1.5f,  1.5f }; // maximum of AABB surrounding the isofunction
    // general settings
    int32_t chebyshevMethod = 0;                    // 0 for sturm sequence based and 1 for qr-based root finding
    bool offsetErrors = false;                      // determine finite step errors, WARNING: very slow
    bool lipschitz = false;                         // evaluate approximate finite lipschitz constant of isofunction, WARNING: slow
    // numerical settings
    bool highPrecisionMode = false;                 // utilize high precision shifting mode as described in the paper
    double lowPrecisionEpsilon = 1.;                // epsilon factor for the second proxy in high precision mode
    double highPrecisionEpsilon = 7.;               // epsilon factor for the initial proxy in high precision mode
    bool polishRoots = true;                        // polish roots using Newton iterations
    int32_t newtonLimit = 16;                       // maximum number of Newton iterations that should be performed
    // misc settings
    int32_t breakPoint = 2;                         // value to shift function by based on approximation errors (with factor 10^breakpoint)

    int32_t intersectionMethod = 1;

    // sphere and segment tracing parameters
    int32_t maximumIterations = -1;
    bool autoLipschitz = true;
    double globalLipschitz = 1.0;
    double segmentAmplification = 1.0;
    double sstThreshold = -5.0;

    // Sherstyuk parameters
    int32_t subIntervalCount = 2;

    // AMP Parameters
    int32_t AMPIntervalCount = 4;
    double tau1 = 1e-3;
    double tau2 = 1e-1;
    double tau3 = 1e-1;
    double AMPepsilon = 1e-12;
    bool useTaylorTest = false;

};
struct isoFunctionStatistics {
    double3 origin;
    double3 direction;
    int32_t rayIndex;
    // hit information
    bool intersectedAABB = false;                   // indicates if the scalar field bounding box was intersected
    bool intersectedSurface = false;                // indicates if the actual isosurface was intersected
    double depth = 0.0;                             // distance value for which the isosurface was intersected
    // chebyshev information
    int32_t degree = 0;                             // degree of base-proxy function
    double minCoeff = 0.0;                          // minimum coefficient of base-proxy function
    double maxCoeff = 0.0;                          // maximum coefficient of base-proxy function
    //std::vector<double> coefficients;               // vector containing all coefficients of base-proxy function
    // proxy information (not always defined)
    bool hasLipschitzInformation = false;
    double fmin, fmax;                              // approximate minimum and maximum function value using base-proxy
    double gmin, gmax;                              // approximate minimum and maximum function derivative value using base-proxy
    // proxy information (always defined)
    double epsilon_a = 0.;                          // degree^2 based approximation error
    double epsilon_u = 0.;                          // cheb::eps based approximation error
    double value = 0.;                              // value that isofunction was shifted by based on approximation errors
    // error information (not always defined)
    bool hasIntegerErrors = false;
    int32_t error = 0;                              // error in finite precision steps on the actual function
    int32_t errorP = 0;                             // error in finite precision steps on the proxy function
    // error information based on function
    double ferror = 0.;                             // result of evaluating actual function at approximate intersection
    double ferrorp = 0.;                            // result of evaluating proxy function at approximate intersection
    // timing information
    double timeForProxy = 0.0;                      // time to construct proxy for the given ray in ms
    double timeForRootFinding = 0.0;                // time to find roots on proxy function
    // surface information
    double3 gradients{ 0.,0.,0. };                  // approximate gradients of scalar field at intersection
    // Chebyshev Surfaces
    int32_t numRoots = 0;                           // approximate number of roots along ray, only accurate without subdivisions & high precision

    // sphere tracing
    int64_t sphereTracingIterations;
    double sphereTracingStep;
    // segment tracing
    int32_t segmentTracingIterations;
    double segmentTracingStep;
    double localLipschitzConstant;
    // LG Surfaces
    int32_t LGdepth;
    int32_t RFdepth;
    double localLConstant;
    double localGConstant;
    // AMP 09
    int32_t AMPIterations;
    double AMPStep;
    double AMPWidth;
    // Sherstyuk
    int32_t SherstyukIterations;
    std::tuple<double, double, double, double, double, double> SherstyukCoefficients;
    // misc information
    bool hasRankInformation = false;
    int3 ranks{ 0,0,0 };                            // if chebyshev gradients were used this contains the degrees of the gradient functions
};


        std::cout << "Reading from " << file << std::endl;
        std::ifstream in;
        in.open(file, std::ios::binary);
        if (!in.is_open()) { std::cerr << "Could not open file!" << std::endl; throw std::invalid_argument{ "File does not exist" }; }

        std::size_t bufferSize = std::max(sizeof(isoFunctionState), sizeof(isoFunctionStatistics));
        std::byte* buffer = new std::byte[bufferSize + 1];
        std::size_t numCoeffs = 0;
        in.read((char*)buffer, sizeof(isoFunctionState));
        //std::cout << sizeof(isoFunctionState) << "\t" << sizeof isoFunctionStatistics << std::endl;

        auto state = *reinterpret_cast<isoFunctionState*>(buffer);
        std::vector<std::vector<double>> coefficients;
        std::vector<isoFunctionStatistics> data;
        data.resize(width * height);
        if (state.intersectionMethod == 1 && state.chebyshevMethod == 1) {
            coefficients.resize(width * height);
        }
        auto pb = progressBar(width * height);
       // std::atomic<int32_t> i = 0;
        for (int32_t x = 0; x < width; ++x) {
            for (int32_t y = 0; y < height; ++y) {
                //pb.update(x, y, 1);
                //pb.update(++i, width);
                in.read((char*)buffer, sizeof(isoFunctionStatistics));
                //std::cout << numCoeffs << std::endl;
                isoFunctionStatistics& pixelStats = *reinterpret_cast<isoFunctionStatistics*>(buffer);
                //std::cout << pixelStats.rayIndex << "\t";
                //std::cout << numCoeffs << std::endl;
                //std::cout << numCoeffs << std::endl;
                auto tempBuff = new std::vector<double>();
                //memcpy((char*) &pixelStats.coefficients, (char*) tempBuff, sizeof(*tempBuff));
                if (state.intersectionMethod == 1 && state.chebyshevMethod == 1) {
                    in.read((char*)&numCoeffs, sizeof(numCoeffs));
                    auto localCoefficients = std::vector<double>();
                    localCoefficients.resize(numCoeffs);
                    if (numCoeffs > 0)
                        in.read((char*)localCoefficients.data(), sizeof(double) * numCoeffs);
                    coefficients[pixelStats.rayIndex] = localCoefficients;
                }
                else {
                    //pixelStats.coefficients = std::vector<double>();
                }
                auto i = (height - y - 1) * width + x;
                data[i] = pixelStats;
            }
        }


    return coefficients;
}


	struct complexState{
		FunctionState r, i;
	};
	complexState deriv(cheb::complex val, const cheb::svec& ak) {
		if (ak.size() == 0)
			return complexState{};
		if (ak.size() == 1)
			return complexState{.r = FunctionState{.f = ak[0]}};

		std::vector<complexState> b(ak.size() + 2);
		auto x = val.real();
		auto y = val.imag();


		for (int32_t k = (int32_t)ak.size() - 1; k >= 0; k -= 1) {	
			auto s = k == 0 ? 1. : 2.;


			b[k].r.f        = ak[k] - b[k+2].r.f        + s * (x * b[k+1].r.f                                            - y * b[k+1].i.f);

			b[k].r.J.dfdx   =       - b[k+2].r.J.dfdx   + s * (x * b[k+1].r.J.dfdx                     + b[k+1].r.f      - y * b[k+1].i.J.dfdx);

			b[k].r.H.d2fdx2 =       - b[k+2].r.H.d2fdx2 + s * (x * b[k+1].r.H.d2fdx2 + b[k+1].r.J.dfdx + b[k+1].r.J.dfdx - y * b[k+1].i.H.d2fdx2);
			b[k].r.H.d2fdxy =       - b[k+2].r.H.d2fdxy + s * (x * b[k+1].r.H.d2fdxy + b[k+1].r.J.dfdy                   - y * b[k+1].i.H.d2fdxy - b[k+1].i.J.dfdx);

			b[k].r.J.dfdy   =       - b[k+2].r.J.dfdy   + s * (x * b[k+1].r.J.dfdy                                       - y * b[k+1].i.J.dfdy - b[k+1].i.f);

			b[k].r.H.d2fdyx =       - b[k+2].r.H.d2fdyx + s * (x * b[k+1].r.H.d2fdyx + b[k+1].r.J.dfdy                   - y * b[k+1].i.H.d2fdyx - b[k+1].i.J.dfdx);
			b[k].r.H.d2fdy2 =       - b[k+2].r.H.d2fdy2 + s * (x * b[k+1].r.H.d2fdy2                                     - y * b[k+1].i.H.d2fdy2 - b[k+1].i.J.dfdy - b[k+1].i.J.dfdy);

			b[k].i.f        =       - b[k+2].i.f        + s * (x * b[k+1].i.f                                            + y * b[k+1].r.f);

			b[k].i.J.dfdx   =       - b[k+2].i.J.dfdx   + s * (x * b[k+1].i.J.dfdx                     + b[k+1].i.f      + y * b[k+1].r.J.dfdx);

			b[k].i.H.d2fdx2 =       - b[k+2].i.H.d2fdx2 + s * (x * b[k+1].i.H.d2fdx2 + b[k+1].i.J.dfdx + b[k+1].i.J.dfdx + y * b[k+1].r.H.d2fdx2);
			b[k].i.H.d2fdxy =       - b[k+2].i.H.d2fdxy + s * (x * b[k+1].i.H.d2fdxy + b[k+1].i.J.dfdy                   + y * b[k+1].r.H.d2fdxy + b[k+1].r.J.dfdx);

			b[k].i.J.dfdy   =       - b[k+2].i.J.dfdy   + s * (x * b[k+1].i.J.dfdy                                       + y * b[k+1].r.J.dfdy   + b[k+1].r.f);

			b[k].i.H.d2fdyx =       - b[k+2].i.H.d2fdyx + s * (x * b[k+1].i.H.d2fdyx + b[k+1].i.J.dfdy                   + y * b[k+1].r.H.d2fdyx + b[k+1].r.J.dfdx);
			b[k].i.H.d2fdy2 =       - b[k+2].i.H.d2fdy2 + s * (x * b[k+1].i.H.d2fdy2                                     + y * b[k+1].r.H.d2fdy2 + b[k+1].r.J.dfdy + b[k+1].r.J.dfdy);
		}

		return b[0];
	}


#include <tools/timer.h>
int main(int argc, char* argv[]) 
{
    // globalCoefficients = chebyshevData("/home/winchenbach/dev/surfaceData/01 - Sphere/Sphere - Chebyshev QR-Decomposition -  - 2022-04-02_14-48-11.log");
    globalCoefficients = chebyshevData("/home/winchenbach/dev/surfaceData/07 - Log SoR/SOR - Chebyshev QR-Decomposition -  - 2022-04-03_13-21-17.log");
    // globalCoefficients = chebyshevData("/home/winchenbach/dev/surfaceData/02 - Cube/Cube - Chebyshev QR-Decomposition -  - 2022-04-03_16-39-04.log");


    std::vector<scalar> coefficients{ 0x1.0000000000000p+0, -0x1.bfb5dbce01c41p-5, 0x1.08a2242d1a359p-1, 0x1.d7d54c363507dp-4, -0x1.5ab0400be2843p-2, -0x1.d1eae3a019217p-4, 0x1.53467431323d9p-3, 0x1.5eca3b803d81fp-4, -0x1.58e81cada5e99p-5, -0x1.e6236ce4c4606p-5, -0x1.c12efc8eac5e8p-6, 0x1.6276d0bf5035ep-5, 0x1.bec289393e0e2p-5, -0x1.1a219e6b86f2cp-5, -0x1.be6f9dca8d40bp-5, 0x1.c86edb95cfa5ep-6, 0x1.4faaec94f4b4cp-5, -0x1.531681d75560ep-6, -0x1.8b071512834abp-6, 0x1.9eec819b561b8p-7, 0x1.4685577420891p-7, -0x1.407ccf235ffb4p-8, -0x1.f826949ba693ap-11, -0x1.4d20028064b9ep-10, -0x1.8f467b5f5dfb1p-9, 0x1.4848b898bdbdcp-8, 0x1.cf8dedbe3c789p-9, -0x1.9b8c05d32b86cp-8, -0x1.2a668ad17ef8fp-9, 0x1.722b800f2f559p-8, 0x1.6163325b83c0dp-11, -0x1.03deb36dc62fdp-8, 0x1.ed102d5493f9cp-12, 0x1.0c0a899298079p-9, -0x1.e0e11a665aa61p-11, -0x1.e1fcde06475d1p-12, 0x1.a17d424cfcff5p-11, -0x1.1420db519d807p-11, -0x1.9620e087d4793p-12, 0x1.debd026d02dd2p-11, -0x1.77e2d9b1cfb88p-15, -0x1.c2a4c616293c8p-11, 0x1.5fa4e33273f1cp-12, 0x1.2f25b30b64a58p-11, -0x1.c728c92780f56p-12, -0x1.09bebd3bcdd60p-12, 0x1.8969dcd9e5390p-12, 0x1.680273ede1100p-20, -0x1.e695886d22a44p-13, 0x1.19d3db31d1c22p-13, 0x1.451d862a51092p-14, -0x1.5830b483f469ep-13, 0x1.47feae7218d71p-15, 0x1.080519441cabcp-13, -0x1.9dbd1095e7c3dp-14, -0x1.040581ab6cdc2p-14, 0x1.bdd7093daa787p-14, 0x1.1edc0b9d64c65p-18, -0x1.51cc6693f858ep-14, 0x1.05395b4abaeebp-15, 0x1.6338e2adf00abp-15, -0x1.613a8e4c9e76cp-15, -0x1.0a23005ae65a4p-17, 0x1.2713dca1a13fcp-15, -0x1.dc2544c9d7998p-17, -0x1.4d0ded3ce256ep-16, 0x1.7b9b2a8f7e71bp-16, 0x1.2590172e562fap-18, -0x1.59f326c541e35p-16, 0x1.a9b4323a1c73fp-18, 0x1.ba5470220865dp-17, -0x1.6d1c5c6a91d44p-17, -0x1.463a52fd827b3p-18, 0x1.5a32cadecb44bp-17, -0x1.80b6e24fbbfefp-20, -0x1.cd39e56e96a2fp-18, 0x1.366773fa23429p-18, 0x1.71cab65c2fec4p-19, -0x1.5205f7fca7e6fp-18, 0x1.14c353fe0d9d4p-21, 0x1.f16b4fcc3e700p-19, -0x1.33f1dd603abf4p-19, -0x1.dbb6cf7610516p-20, 0x1.64084d89122d9p-19, 0x1.6dd3f47632fabp-24, -0x1.150322e20639cp-19, 0x1.f5c7d0c7726b0p-21, 0x1.26036598cdf96p-20, -0x1.4ec6fff1a3e10p-20, -0x1.a61dad27fe1dap-23, 0x1.190d666bced26p-20, -0x1.a1469b963d869p-22, -0x1.45a66f0c91cf0p-21, 0x1.48ad8898fca18p-21, 0x1.5796bcdd57b4cp-23, -0x1.2a75cb7178f83p-21, 0x1.489cc531d2b5cp-23, 0x1.7af0d6daebba4p-22, -0x1.38b650ef06663p-22, -0x1.0e96f5eaffd60p-23, 0x1.30eb66ea77323p-22, -0x1.86dd2765c282fp-25, -0x1.9eb96c1b94f0dp-23, 0x1.1b87f08109b25p-23, 0x1.5ab30fd3869f0p-24, -0x1.316dec25e342cp-23, 0x1.59ac252061764p-27, 0x1.c0f26c684ecb2p-24, -0x1.037615fb53c96p-24, -0x1.aff1985341829p-25, 0x1.3385568a55678p-24, 0x1.975029396f7d1p-29, -0x1.e3f215e378686p-25, 0x1.bada768c2a27dp-26, 0x1.03418d1154959p-25, -0x1.2acce3e83b422p-25, -0x1.812e5674328bep-28, 0x1.f9930aa1acfeep-26, -0x1.6a0b052bdd143p-27, -0x1.289a48b8d3084p-26, 0x1.207e243827d7fp-26, 0x1.4b580bdc3afa5p-28, -0x1.075a94e0fe19ap-26, 0x1.113de10ebf786p-28, 0x1.5006ec8f4406dp-27, -0x1.10b92ed181238p-27, -0x1.e4cdf13420658p-29, 0x1.0d6b308d1b188p-27, -0x1.50acc4addcea3p-30, -0x1.723cbcdc933dcp-28, 0x1.f48281f359180p-29, 0x1.3bd4524d638a6p-29, -0x1.0f3992e52404dp-28, 0x1.e66b0c53c1961p-33, 0x1.910480e16b5f9p-29, -0x1.c152ed8ced550p-30, -0x1.863b4865a9378p-30, 0x1.0e2d141dd36aap-29, 0x1.cb1835183e3ccp-34, -0x1.aca6f786691b9p-30 };
    coefficients = globalCoefficients[0];

     globalFunction = cheb::Function(std::vector<cheb::IntervalFunction>{cheb::IntervalFunction(coefficients)});

    // globalFunction = cheb::Function([](cheb::scalar x){ return std::cos(x * 13.) * sin(x * 25.);},cheb::Domain{-1.,1.});
    // globalFunction = cheb::Function(std::vector<cheb::IntervalFunction>{cheb::IntervalFunction(coefficients)});
    // globalFunction = cheb::Function([](cheb::scalar x) { return std::tan(x); }, cheb::Domain{-1.,1.});
    // globalFunction = cheb::Function([](cheb::scalar x) { return x * x * x - 1.; }, cheb::Domain{ -1.,1. });
    // globalFunction = cheb::Function([](cheb::scalar x) { return x * x * x - 2. * x + 2.; }, cheb::Domain{ -1.,1. });
    // globalFunction = cheb::Function([](cheb::scalar x) { return x * x * x * x * x - 1; }, cheb::Domain{ -1.,1. });

    // globalFunction = cheb::Function([](cheb::scalar x) { return std::exp(x) * std::sin(x); }, cheb::Domain{-1.,1.});



    globalFunctionFirstDerivative = globalFunction.derivative();
    globalFunctionSecondDerivative = globalFunctionFirstDerivative.derivative();

    std::cout << evalFunction(cheb::complex(1., 0.)) << " :8 " << globalFunction(1.) << std::endl;

    for (auto c : globalFunction.funs[0].coeffs())
        std::cout << c << " ";
    std::cout << std::endl;

    auto eval = [](cheb::complex location){
        auto fn = [](scalar x, scalar y){
                cheb::complex location(x,y);
            auto [fr, fi] = deriv(location, globalFunction.funs[0].coeffs());

            auto f = fr.f * fr.f + fi.f * fi.f;

            auto dfdx = 2. * fr.f * fr.J.dfdx + 2. * fi.f * fi.J.dfdx;
            auto dfdy = 2. * fr.f * fr.J.dfdy + 2. * fi.f * fi.J.dfdy;

            auto d2fdx2 = 2. * (fr.J.dfdx * fr.J.dfdx + fr.f * fr.H.d2fdx2 + fi.J.dfdx * fi.J.dfdx + fi.f * fi.H.d2fdx2);
            auto d2fdxy = 2. * (fr.J.dfdy * fr.J.dfdx + fr.f * fr.H.d2fdxy + fi.J.dfdy * fi.J.dfdx + fi.f * fi.H.d2fdxy);
            auto d2fdyx = 2. * (fr.J.dfdy * fr.J.dfdx + fr.f * fr.H.d2fdyx + fi.J.dfdy * fi.J.dfdx + fi.f * fi.H.d2fdyx);
            auto d2fdy2 = 2. * (fr.J.dfdy * fr.J.dfdy + fr.f * fr.H.d2fdy2 + fi.J.dfdy * fi.J.dfdy + fi.f * fi.H.d2fdy2);

            Jacobian J{dfdx,dfdy};
            Hessian H{d2fdx2, d2fdxy, d2fdyx, d2fdy2};
	        return FunctionState{f,J,H};
        };
        auto [fr, fi] = deriv(location, globalFunction.funs[0].coeffs());

        auto x = location.real();
        auto y = location.imag();
        auto [f,J,H] = fn(x,y);
        auto h = 1e-6;
        auto [fx, Jx, Hx] = fn(x + h, y);
        auto [fy, Jy, Hy] = fn(x, y + h);


        std::cout << "f( " << location.real() << " + " << location.imag() << "i ) = " << fr.f << " + " << fi.f << "i -> " << f << std::endl;

        std::cout << J.dfdx << " : " << J.dfdy << " | " << (fx - f) / h << " : " << (fy - f) / h << std::endl;

        std::cout << H.d2fdx2 << " | " << (Jx.dfdx - J.dfdx) / h << std::endl;
        std::cout << H.d2fdxy << " | " << (Jx.dfdy - J.dfdy) / h << std::endl;
        std::cout << H.d2fdyx << " | " << (Jy.dfdx - J.dfdx) / h << std::endl;
        std::cout << H.d2fdy2 << " | " << (Jy.dfdy - J.dfdy) / h << std::endl;

        

        auto det = H.d2fdx2 * H.d2fdy2 - H.d2fdxy * H.d2fdyx;
        Hessian H_1{ H.d2fdy2 / det, -H.d2fdxy / det, -H.d2fdyx / det, H.d2fdx2/det};

        auto prodx = H_1.d2fdx2 * J.dfdx + H_1.d2fdxy * J.dfdy;
        auto prody = H_1.d2fdyx * J.dfdx + H_1.d2fdy2 * J.dfdy;

        auto dx = cheb::complex(prodx, prody);
        std::cout << det << " : " << H_1.d2fdx2 << " " << H_1.d2fdxy << " " << H_1.d2fdyx << " " << H_1.d2fdy2 << std::endl;

        std::cout << "[ " << prodx << " : " << prody << " ] -- [ " << -J.dfdx << " : " << -J.dfdy << " ]" << std::endl;


        std::cout<< std::endl;
    };

    eval(cheb::complex(-1, -1));
    eval(cheb::complex(-1, 1));
    eval(cheb::complex(1, -1));
    eval(cheb::complex(1, 1));
    eval(cheb::complex(0, 0));
    //system("pause");
   // int x = getchar();


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