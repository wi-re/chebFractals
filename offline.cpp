#include "gui/glui.h"
#include <sstream>
#include <thread>
#include <iostream> 
#include <cheb/cheb.h>
#include <iomanip>
#include <vector>
#include <fstream>
#include <filesystem>
#include <tools/stb_image_write.h>
#include <tools/stb_image.h>
#include <mutex>
#include <omp.h>
#include <cfloat>

struct progressBar {
    std::chrono::high_resolution_clock::time_point start, lastTime;
    int32_t frames, lastFrame = 0;

    progressBar(int32_t targetFrames) : frames(targetFrames) { start = std::chrono::high_resolution_clock::now(); };

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
                auto totalTime = ((float)std::chrono::duration_cast<std::chrono::microseconds>(now - start).count()) / 1000.f / 1000.f;
                std::cout << std::fixed << std::setprecision(0) << " " << std::setw(3) << totalTime << "/" << std::setw(3) << (totalTime / progress) << "s  ";
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
                auto totalTime = ((float)std::chrono::duration_cast<std::chrono::microseconds>(now - start).count()) / 1000.f / 1000.f;
                std::cout << std::fixed << std::setprecision(0) << " " << std::setw(3) << totalTime << "/" << std::setw(3) << (totalTime / progress) << "s  ";
            }
        }
        std::cout << "\r";
        std::cout.flush();
        std::cout.copyfmt(cout_state);
    }
    void end() { std::cout << std::endl; }
};


auto getColorMap() {
    auto colorMap = ParameterManager::instance().get<std::string>("colorMap.colorMap");
    int32_t image_width = 1024;
    int32_t image_height = 1024;
    v4* img = new v4[1024];
    v4* color_map = nullptr;
    for (int32_t it = 0; it < 1024; ++it)
        img[it] = v4{ (float)it / (float)1024 * 255.f, (float)it / (float)1024 * 255.f, (float)it / (float)1024 * 255.f, 255.f };
    // std::string file_name = std::string("cfg/") + colorMap + ".png";
    try {
        std::string file_name = resolveFile(std::string("cfg/") + ParameterManager::instance().get<std::string>("colorMap.colorMap") + ".png").string();

        if (std::filesystem::exists(file_name)) {
            std::cout << "Loading " << file_name << std::endl;
            unsigned char* image_data = stbi_load(file_name.c_str(), &image_width, &image_height, NULL, 4);
            delete[] img;
            img = new v4[image_width];
            for (int32_t it = 0; it < image_width; ++it) {
                img[it] = v4{ (float)image_data[it * 4 + 0 + image_width * (image_height / 2)], (float)image_data[it * 4 + 1 + image_width * (image_height / 2)],
                             (float)image_data[it * 4 + 2 + image_width * (image_height / 2)], 255.f };
                // std::cout << it << ": [ " << img[it].x << " " << img[it].y << " "
                // << img[it].z <<
                //     " ]" << std::endl;
            }
            // img = QImage(QString::fromStdString(file_name));
            // img.load(QString(file_name.c_str()));
            // std::cout << image_width << " : " << image_height << std::endl;
        }
    }
    catch (...) {
    }
    color_map = (v4*)realloc(color_map, sizeof(v4) * (image_width));
    for (int32_t it = 0; it < image_width; ++it) {
        color_map[it] = v4{ (float)(img[it].x) / 256.f, (float)(img[it].y) / 256.f, (float)(img[it].z) / 256.f, 1.f };
        // std::cout << color_map[it].x() << " : " << color_map[it].y() << " : "
        // << color_map[it].z() << std::endl; if(it == img.width()-1)
        //	color_map[it + 1] = QVector4D{ (float)(col.red()) / 256.f,
        //(float)(col.green()) / 256.f, 	(float)(col.blue()) / 256.f, 1.f };
    }
    auto color_map_elements = image_width;

    return std::make_tuple(colorMap, color_map, color_map_elements);
}

auto mapToCmap(double val, double min, double max, v4* cmap, int32_t cmapElements) {
    auto xrel = std::clamp((val - min) / (max - min), 0., 1.);

    auto dx = 1. / (double)(cmapElements - 1);

    auto left = std::clamp((int32_t)::floor(xrel / dx),0,cmapElements - 1);
    auto right = std::clamp((int32_t)::ceil(xrel / dx), 0, cmapElements - 1);

    auto alpha = 1.f - (float) ((xrel - ((double)left) * dx) / dx);

    auto vLeft = cmap[left];
    auto vRight = cmap[right];

    return v4{
        vLeft.x * alpha + vRight.x * (1.f - alpha),
        vLeft.y * alpha + vRight.y * (1.f - alpha),
        vLeft.z * alpha + vRight.z * (1.f - alpha),
        vLeft.w * alpha + vRight.w * (1.f - alpha) };
}



int main(int argc, char* argv[]) 
{
    int32_t batchSize = 1024;

    auto& gui = GUI::instance();
    gui.render_lock.lock();
    gui.initParameters(argc, argv);
    gui.initGVDB();
    gui.initSimulation();   

    ParameterManager::instance().load("render.yaml");

    const std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
    const std::time_t t_c = std::chrono::system_clock::to_time_t(now);
    auto time = std::time(nullptr);
    std::stringstream ss;
    char buf[256];
    ss << std::put_time(std::localtime(&t_c), "%F_%T") << "_" << ParameterManager::instance().get<std::string>("field.method")
        << ".png"; // ISO 8601 without timezone information.
    auto s = ss.str();
    std::replace(s.begin(), s.end(), ':', '-');


    auto [colorMap, cmap, cmapElements] = getColorMap();

    std::ifstream config("render.config");
    static auto& xmin = ParameterManager::instance().get<scalar>("domain.xmin");
    static auto& xmax = ParameterManager::instance().get<scalar>("domain.xmax");
    static auto& ymin = ParameterManager::instance().get<scalar>("domain.ymin");
    static auto& ymax = ParameterManager::instance().get<scalar>("domain.ymax");
    std::string t;

    config >> t >> std::hexfloat >> xmin;
    std::cout << t << " " << xmin << std::endl;
    config >> t >> std::hexfloat >> ymin;
    std::cout << t << " " << ymin << std::endl;
    config >> t >> std::hexfloat >> xmax;
    std::cout << t << " " << xmax << std::endl;
    config >> t >> std::hexfloat >> ymax;
    std::cout << t << " " << ymax << std::endl;
   
    int32_t baseWidth, baseHeight, multiplier, coeffs;

    config >> t >> baseWidth;
    std::cout << t << " " << baseWidth << std::endl;
    config >> t >> baseHeight;
    std::cout << t << " " << baseHeight << std::endl;
    config >> t >> multiplier;
    std::cout << t << " " << multiplier << std::endl;
    config >> t >> coeffs;
    std::cout << t << " " << coeffs << std::endl;

    std::vector<double> coefficients;

    while (config) {
        double c;
        config >> c;
        coefficients.push_back(c);
    }
    config.close();
    std::cout << coefficients.size() << " " << coeffs << std::endl;
    globalFunction = cheb::Function(std::vector<cheb::IntervalFunction>{cheb::IntervalFunction(coefficients)});

    int32_t xbatches = baseWidth * multiplier / batchSize + ((baseWidth * multiplier) % batchSize == 0 ? 0 : 1);
    int32_t ybatches = baseHeight * multiplier / batchSize +( (baseHeight * multiplier) % batchSize == 0 ? 0 : 1);

    int32_t xLim = baseWidth * multiplier;
    int32_t yLim = baseHeight * multiplier;

    std::cout << "Rendering to " << xbatches << " x " << ybatches << " @ " << batchSize << " x " << batchSize << "px";

    std::filesystem::create_directories("renderOutput");


    static int32_t* buffer = new int32_t[batchSize * batchSize];
    char* buffC = (char*)buffer;
    static char* bufferC = new char[batchSize * batchSize * 4];
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist(0.f, 1.f);
    auto stringMethod = ParameterManager::instance().get<std::string>("field.method");
    bool cycles = ParameterManager::instance().get<bool>("field.cycles");
    optimizationMethod method = getMethod(stringMethod);

    std::cout << "max threads: " << omp_get_max_threads() << std::endl << std::endl;
    for (int32_t xi = 0; xi < xbatches; ++xi) {
        for (int32_t yi = 0; yi < ybatches; ++yi) {
            std::cout << "Processing batch " << xi << " x " << yi << std::endl;
            std::stringstream s;
            s << "renderOutput" << "/" << std::setw(4) << std::setfill('0') << yi << "_" << std::setw(4) << std::setfill('0') << xi << ".png";

            if (std::filesystem::exists(s.str())) {
                std::cout << "File already exists" << std::endl;
            }

            int32_t pxmin = xi * batchSize;
            int32_t pxmax = xi * batchSize + batchSize;

            int32_t pymin = yi * batchSize;
            int32_t pymax = yi * batchSize + batchSize;

            std::cout << "Pixel range: [ " << pxmin << " : " << pymin << " ] -> [ " << pxmax << " : " << pymax << " ]" << std::endl;

            progressBar pbl{ batchSize * batchSize };
            std::atomic<int32_t> iterationCounter = 0;

            #pragma omp parallel for schedule(dynamic, 32) num_threads(32)
            for (int32_t j = 0; j < batchSize; ++j) {
                int32_t actualy = pymin + j;
                for (int32_t i = 0; i < batchSize; ++i) {
                    int32_t actualx = pxmin + i;
                    if (actualy >= yLim || actualx >= xLim) {
                        buffer[j * batchSize + i] = 0;
                    }
                    else {
                        double xRel = ((double)actualx) / ((double)xLim);
                        double yRel = ((double)actualy) / ((double)yLim);
                        double xPos = xmin + xRel * (xmax - xmin);
                        double yPos = ymin + yRel * (ymax - ymin);

                        cheb::complex pos = cheb::complex(xPos, yPos);

                        realOffset = ParameterManager::instance().get<scalar>("field.offset");
                        complexOffset = ParameterManager::instance().get<scalar>("field.coffset");

                        auto [state, location, positions, values, steps] = optimize(pos, method);

                        auto angle = std::arg(location);
                        auto rad = std::abs(location);
                        v4 mapped = mapToCmap(rad > 1e-3 ? angle : std::numbers::pi, -std::numbers::pi, std::numbers::pi, cmap, cmapElements);

                        buffC[(j * batchSize + i) * 4 + 0] = mapped.x * 255.f;
                        buffC[(j * batchSize + i) * 4 + 1] = mapped.y * 255.f;
                        buffC[(j * batchSize + i) * 4 + 2] = mapped.z * 255.f;
                        buffC[(j * batchSize + i) * 4 + 3] = 255;
                        
                        if (state == optimizationState::diverged) {
                            buffC[(j * batchSize + i) * 4 + 0] = 255;
                            buffC[(j * batchSize + i) * 4 + 1] = 0;
                            buffC[(j * batchSize + i) * 4 + 2] = 0;
                            buffC[(j * batchSize + i) * 4 + 3] = 255;
                        }
                        else if(cycles && state == optimizationState::unstable) {
                            buffC[(j * batchSize + i) * 4 + 0] = 0;
                            buffC[(j * batchSize + i) * 4 + 1] = 255;
                            buffC[(j * batchSize + i) * 4 + 2] = 0;
                            buffC[(j * batchSize + i) * 4 + 3] = 255;
                        }
                    }
                    int32_t iter = ++iterationCounter;
                    if (iter % (4 * 8192) == 0) {
                    #pragma omp critical
                        { pbl.update(iter, 1); }
                    }
                }
            }
            for (int32_t i = 0; i < batchSize; i++) {
                for (int32_t j = 0; j < batchSize; j++) {
                    bufferC[(i + batchSize * (batchSize - j - 1)) * 4 + 0] = buffC[(i + batchSize * j) * 4 + 0];
                    bufferC[(i + batchSize * (batchSize - j - 1)) * 4 + 1] = buffC[(i + batchSize * j) * 4 + 1];
                    bufferC[(i + batchSize * (batchSize - j - 1)) * 4 + 2] = buffC[(i + batchSize * j) * 4 + 2];
                    bufferC[(i + batchSize * (batchSize - j - 1)) * 4 + 3] = buffC[(i + batchSize * j) * 4 + 3];
                }
            }


            if (!stbi_write_png(s.str().c_str(), batchSize, batchSize, 4, bufferC, batchSize * 4)) {
                std::cerr << "ERROR: could not write image to " << s.str() << std::endl;
            }


        }
    }

    return 0;
}