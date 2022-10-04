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

#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"
#include "imgui/implot.h"

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
GLFWwindow* window;
int32_t* imageBuffer;

int64_t xLim = 0;
int64_t yLim = 0;

double data_xmin = 0.;
double data_xmax = 0.;
double data_ymin = 0.;
double data_ymax = 0.;

void OSD() {

    auto [minDomain, maxDomain] = getDomain();
    double xp, yp;
    glfwGetCursorPos(window, &xp, &yp);
    vec m_cursorPosition(xp,yp);

    auto x = std::clamp(m_cursorPosition.x(), 0.0, (scalar)screenWidth) / (scalar)screenWidth;
    auto y = std::clamp(m_cursorPosition.y(), 0.0, (scalar)screenHeight) / (scalar)screenHeight;
    if (y > 1.0)
        y -= 1.0;
    y = 1.0 - y;
    x *= maxDomain.x() - minDomain.x();
    y *= maxDomain.y() - minDomain.y();
    x += minDomain.x();
    y += minDomain.y();

    vec position(x, y);

    const float DISTANCE = 10.0f;
    static int corner = 3;
    ImGuiIO& io = ImGui::GetIO();
    if (corner != -1) {
        ImGuiViewport* viewport = ImGui::GetMainViewport();
        ImVec2 work_area_pos = viewport->WorkPos; // Instead of using viewport->Pos we use GetWorkPos() to avoid menu bars, if any!
        ImVec2 work_area_size = viewport->WorkSize;
        ImVec2 window_pos = ImVec2((corner & 1) ? (work_area_pos.x + work_area_size.x - DISTANCE) : (work_area_pos.x + DISTANCE),
            (corner & 2) ? (work_area_pos.y + work_area_size.y - DISTANCE) : (work_area_pos.y + DISTANCE));
        ImVec2 window_pos_pivot = ImVec2((corner & 1) ? 1.0f : 0.0f, (corner & 2) ? 1.0f : 0.0f);
        ImGui::SetNextWindowPos(window_pos, ImGuiCond_Always, window_pos_pivot);
        ImGui::SetNextWindowViewport(viewport->ID);
    }
    ImGui::SetNextWindowBgAlpha(0.35f); // Transparent background
    bool m_showInfo = true;
    if (ImGui::Begin("Particle Info", &m_showInfo,
        (corner != -1 ? ImGuiWindowFlags_NoMove : 0) | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
        ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav)) {
        auto addScalar = [&](auto name, scalar value, auto description) {
            float vcp = value;
            auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
            ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
            ImGui::DragFloat(name, &vcp, 0, vcp, vcp, "%g");
            ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("%s", description);
            return;
        };
        struct float2 {
            float x, y;
        };
        auto addScalar2 = [&](auto name, vec value, auto description) {
            float2 vcp(value.x(), value.y());
            auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
            ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
            ImGui::DragFloat2(name, &vcp.x, 0, 0, 0, "%g");
            ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("%s", description);
            return;
        };
        auto addInteger = [&](auto name, int32_t value, auto description) {
            int32_t vcp = value;
            auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
            ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
            ImGui::DragInt(name, &vcp, 0, vcp, vcp);
            ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("%s", description);
            return;
        };

        addScalar2("Cursor        ", position, "");
        addScalar2("Visible Min ", minDomain, "");
        addScalar2("Visible Max ", maxDomain, "");

        vec bottomleft(minDomain.x(), minDomain.y());
        vec topright(maxDomain.x(), maxDomain.y());

        vec data_bottomleft(data_xmin, data_ymin);
        vec data_topright(data_xmax, data_ymax);

        //std::cout << "Visible Range: " << " [ " << bottomleft.x() << " " << bottomleft.y() << " ] -> [ " << topright.x() << " " << topright.y() << " ]\n";
        //std::cout << "Data Range: " << " [ " << data_bottomleft.x() << " " << data_bottomleft.y() << " ] -> [ " << data_topright.x() << " " << data_topright.y() << " ]\n";

        vec rel_bl = (bottomleft - data_bottomleft).cwiseQuotient(data_topright - data_bottomleft);
        vec rel_tr = (topright - data_bottomleft).cwiseQuotient(data_topright - data_bottomleft);

        //std::cout << "Relative: " << " [ " << rel_bl.x() << " " << rel_bl.y() << " ] -> [ " << rel_tr.x() << " " << rel_tr.y() << " ]\n";

        int32_t min_x = (int32_t)(std::clamp(rel_bl.x(), 0., 1.) * (double)xLim);
        int32_t max_x = (int32_t)(std::clamp(rel_tr.x(), 0., 1.) * (double)xLim);
        int32_t min_y = (int32_t)(std::clamp(rel_bl.y(), 0., 1.) * (double)yLim);
        int32_t max_y = (int32_t)(std::clamp(rel_tr.y(), 0., 1.) * (double)yLim);

        //std::cout << "Pixels: " << " [ " << min_x << " " << min_y << " ] -> [ " << max_x << " " << max_y << " ]\n";

        vec visiblePixels(max_x - min_x,max_y - min_y);
        addScalar2("Pixels", visiblePixels, "");


        {
            auto ny = ParameterManager::instance().get<int32_t>("field.ny");
            auto nx = ParameterManager::instance().get<int32_t>("field.nx");
            auto x = std::clamp(m_cursorPosition.x(), 0.0, (scalar)screenWidth) / (scalar)screenWidth;
            auto y = 1. - std::clamp(m_cursorPosition.y(), 0.0, (scalar)screenHeight) / (scalar)screenHeight;

            int32_t xi = std::clamp((int32_t)::floor(x * dataWidth), 0, nx - 1);
            int32_t yi = std::clamp((int32_t)::floor(y * dataHeight), 0, ny - 1);

            addScalar2("Cursor (screen)", vec(x, y), "");
            addScalar2("Cursor (int)", vec(xi, yi), "");
            //scalar value = scalarFieldData[xi + yi * dataWidth];
            //addScalar("Scalar Field", value, "");
        }
    }
    ImGui::End();
}


void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    static bool emit = false;
    static auto& xmin = ParameterManager::instance().get<scalar>("domain.xmin");
    static auto& xmax = ParameterManager::instance().get<scalar>("domain.xmax");
    static auto& ymin = ParameterManager::instance().get<scalar>("domain.ymin");
    static auto& ymax = ParameterManager::instance().get<scalar>("domain.ymax");
    if (action != GLFW_PRESS) {
        auto scaling = 1.05;
        if (mods & GLFW_MOD_SHIFT)
            scaling = 1.01;
        if (mods & GLFW_MOD_CONTROL)
            scaling = 1.5;
        auto dx = (xmax - xmin);
        auto dy = ymax - ymin;
        auto xcenter = (xmax + xmin) / 2.;
        auto ycenter = (ymax + ymin) / 2.;

        switch (key) {
        case GLFW_KEY_LEFT:
            xmin -= dx * (scaling - 1.);
            xmax -= dx * (scaling - 1.);
            break;
        case GLFW_KEY_RIGHT:
            xmin += dx * (scaling - 1.);
            xmax += dx * (scaling - 1.);
            break;
        case GLFW_KEY_UP:
            ymin += dy * (scaling - 1.);
            ymax += dy * (scaling - 1.);
            break;
        case GLFW_KEY_DOWN:
            ymin -= dy * (scaling - 1.);
            ymax -= dy * (scaling - 1.);
            break;
        case GLFW_KEY_APOSTROPHE:
            xmin = xcenter - dx / 2. * scaling;
            xmax = xcenter + dx / 2. * scaling;
            ymin = ycenter - dy / 2. * scaling;
            ymax = ycenter + dy / 2. * scaling;
            break;
        case GLFW_KEY_RIGHT_BRACKET:
            xmin = xcenter - dx / 2. / scaling;
            xmax = xcenter + dx / 2. / scaling;
            ymin = ycenter - dy / 2. / scaling;
            ymax = ycenter + dy / 2. / scaling;
            break;
        }

        return;
    }
    switch (key) {
    case GLFW_KEY_O: {
        auto [minDomain, maxDomain] = getDomain();

        vec bottomleft(minDomain.x(), minDomain.y());
        vec topright(maxDomain.x(), maxDomain.y());

        vec data_bottomleft(data_xmin, data_ymin);
        vec data_topright(data_xmax, data_ymax);

        //std::cout << "Visible Range: " << " [ " << bottomleft.x() << " " << bottomleft.y() << " ] -> [ " << topright.x() << " " << topright.y() << " ]\n";
        //std::cout << "Data Range: " << " [ " << data_bottomleft.x() << " " << data_bottomleft.y() << " ] -> [ " << data_topright.x() << " " << data_topright.y() << " ]\n";

        vec rel_bl = (bottomleft - data_bottomleft).cwiseQuotient(data_topright - data_bottomleft);
        vec rel_tr = (topright - data_bottomleft).cwiseQuotient(data_topright - data_bottomleft);

        //std::cout << "Relative: " << " [ " << rel_bl.x() << " " << rel_bl.y() << " ] -> [ " << rel_tr.x() << " " << rel_tr.y() << " ]\n";

        int32_t min_x = (int32_t)(std::clamp(rel_bl.x(), 0., 1.) * (double)xLim);
        int32_t max_x = (int32_t)(std::clamp(rel_tr.x(), 0., 1.) * (double)xLim);
        int32_t min_y = (int32_t)(std::clamp(rel_bl.y(), 0., 1.) * (double)yLim);
        int32_t max_y = (int32_t)(std::clamp(rel_tr.y(), 0., 1.) * (double)yLim);

        std::size_t width = max_x - min_x;
        std::size_t height = max_y - min_y;

        int32_t* output_buffer = new int32_t[width * height];
        std::cout << "Processing Data" << std::endl;
#pragma omp parallel for
        for(int64_t j = 0; j < height; ++j){
            for (std::size_t i = 0; i < width; ++i) {
                output_buffer[(height - j - 1ull) * width + i] = imageBuffer[min_x + i + (min_y + j) * xLim];
            }
        }
        std::cout << "Processing Done" << std::endl;

        std::string s = "rendered.png";
        std::cout << "Writing Data" << std::endl;

        if (!stbi_write_png(s.c_str(), width, height, 4, (char*)output_buffer, width * 4)) {
            std::cerr << "ERROR: could not write image to " << s << std::endl;
        }
        std::cout << "Done." << std::endl;

        delete[] output_buffer;

    }
        break;
    }
}
void sizeCallback(GLFWwindow* window, int width, int height) {
    int32_t oldWidth = screenWidth;
    int32_t oldHeight = screenHeight;

    screenWidth = width;
    screenHeight = height;

    static auto& xmin = ParameterManager::instance().get<scalar>("domain.xmin");
    static auto& xmax = ParameterManager::instance().get<scalar>("domain.xmax");
    static auto& ymin = ParameterManager::instance().get<scalar>("domain.ymin");
    static auto& ymax = ParameterManager::instance().get<scalar>("domain.ymax");

    auto xgrowth = 1. + ((((double)screenWidth) / ((double)oldWidth) - 1.) / 2.);
    auto ygrowth = 1. + ((((double)screenHeight) / ((double)oldHeight) - 1.) / 2.);

    auto xcenter = (xmin + xmax) / 2.;
    auto ycenter = (ymin + ymax) / 2.;
    auto dx = (xmax - xmin) / 2.;
    auto dy = (ymax - ymin) / 2.;


    std::cout << "[ " << xmin << " : " << xmax << " ] growth: " << xgrowth << " center: " << xcenter << " d: " << dx << std::endl;
    std::cout << "[ " << ymin << " : " << ymax << " ] growth: " << ygrowth << " center: " << ycenter << " d: " << dy << std::endl;

    xmin = xcenter - dx * xgrowth;
    xmax = xcenter + dx * xgrowth;
    ymin = ycenter - dy * ygrowth;
    ymax = ycenter + dy * ygrowth;
}


static void glfw_error_callback(int error, const char* description) { fprintf(stderr, "Glfw Error %d: %s\n", error, description); }
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

    data_xmin = xmin;
    data_xmax = xmax;
    data_ymin = ymin;
    data_ymax = ymax;

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

    xLim = baseWidth * multiplier;
    yLim = baseHeight * multiplier;

    std::cout << "Rendering to " << xbatches << " x " << ybatches << " @ " << batchSize << " x " << batchSize << "px";

    std::filesystem::create_directories("renderOutput");

    screenWidth = baseWidth;
    screenHeight = baseHeight;
    static int32_t* buffer = new int32_t[batchSize * batchSize];
    char* buffC = (char*)buffer;
    static char* bufferC = new char[batchSize * batchSize * 3];
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist(0.f, 1.f);
    auto stringMethod = ParameterManager::instance().get<std::string>("field.method");
    bool cycles = ParameterManager::instance().get<bool>("field.cycles");
    optimizationMethod method = getMethod(stringMethod);

    std::cout << "Allocating Memory with " << xLim * yLim * 4 << "Bytes" << std::endl;
    imageBuffer = new int32_t[xLim * yLim];
    std::cout << "Reading Data" << std::endl;
    stbi_set_flip_vertically_on_load(1);
#pragma omp parallel for
    for (int64_t xi = 0; xi < xbatches; ++xi) {
        for (int64_t yi = 0; yi < ybatches; ++yi) {
            printf("Processing %d x %d \n",xi,yi);
            //std::cout << "Processing batch " << xi << " x " << yi << std::endl;
            std::stringstream s;
            s << "renderOutput" << "/" << std::setw(4) << std::setfill('0') << yi << "_" << std::setw(4) << std::setfill('0') << xi << ".png";

            if (std::filesystem::exists(s.str())) {
                int32_t batchWidth, batchHeight;
                unsigned char* image_data = stbi_load(s.str().c_str(), &batchWidth, &batchHeight, NULL, 4);
                int32_t* img_data4 = (int32_t*)image_data;
//#pragma omp parallel for
                for (int64_t j = 0; j < batchHeight; ++j) {
                    int64_t actual_y = yi * batchHeight + j;
                    for (int64_t i = 0; i < batchWidth; ++i) {
                        int64_t actual_x = xi * batchWidth + i;
                        if (actual_x >= xLim || actual_y >= yLim)
                            break;
                        imageBuffer[xLim * actual_y + actual_x] = img_data4[batchWidth * j + i];
                    }
                    if (actual_y >= yLim)
                        break;
                }
                stbi_image_free(image_data);
                //free(image_data);
            }
        }
    }


    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
        throw std::runtime_error("Could not setup glfw context");
    const char* glsl_version = "#version 150";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);

    window = glfwCreateWindow(screenWidth, screenHeight, "OmniFlow", NULL, NULL);
    if (window == NULL)
        throw std::runtime_error("Could not setup window");
    glfwSetWindowPos(window, 40, 50);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    bool err = gladLoadGL() == 0;
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
    if (err) {
        fprintf(stderr, "Failed to initialize OpenGL loader!\n");
        throw std::runtime_error("Failed to initialize OpenGL loader");
    }
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;
    ImGui::StyleColorsLight();
    ImGuiStyle& style = ImGui::GetStyle();
    if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable) {
        style.WindowRounding = 0.0f;
        style.Colors[ImGuiCol_WindowBg].w = 1.0f;
    }
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);
    ImGui::PushStyleVar(ImGuiStyleVar_ChildRounding, 5.0f);
    io.Fonts->AddFontFromFileTTF("Cascadia.ttf", 16.0f);
    io.Fonts->AddFontFromFileTTF("NotoMono-Regular.ttf", 16.0f);
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    glfwSetWindowSizeCallback(window, [](GLFWwindow* window, int width, int height) {sizeCallback(window, width, height); });

    glfwSetKeyCallback(window, [](GLFWwindow* window, int key, int scancode, int action, int mods) {keyCallback(window, key, scancode, action, mods); });
    //glfwSetKeyCallback(window, [](GLFWwindow* window, int key, int scancode, int action, int mods) { GUI::instance().keyCallback(window, key, scancode, action, mods); });
    //glfwSetCursorPosCallback(window, [](GLFWwindow* window, double xpos, double ypos) { GUI::instance().cursorPositionCallback(window, xpos, ypos); });
    //glfwSetMouseButtonCallback(window, [](GLFWwindow* window, int button, int action, int mods) { GUI::instance().mouseButtonCallback(window, button, action, mods); });
    //   glfwSetScrollCallback(window, [](GLFWwindow *window, double xpos, double ypos) { GUI::instance().scrollCallback(window, xpos, ypos); });
    glfwSwapInterval(0);

    //ImGuiIO& io = ImGui::GetIO();
    //(void)io;


    auto vtxShader = R"(#version 150

// Input vertex data, different for all executions of this shader.
in vec3 vertexPosition_modelspace;

// Output data ; will be interpolated for each fragment.
out vec2 UV;

void main(){
	gl_Position =  vec4(vertexPosition_modelspace,1);
	UV = (vertexPosition_modelspace.xy+vec2(1,1))/2.0;
})";
    auto frgShader = R"(#version 150

in vec2 UV;

out vec3 color;

uniform sampler2D renderedTexture;

const float PI = 3.1415926535897932384626433832795;

vec3 hsl2rgb( vec3 c )
{
    vec3 rgb = clamp( abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),6.0)-3.0)-1.0, 0.0, 1.0 );

    return c.z + c.y * (rgb-0.5)*(1.0-abs(2.0*c.z-1.0));
}

void main(){
vec4 v = texture( renderedTexture, UV );
color = vec3(v.xyz);
})";
    GLuint quad_VertexArrayID;
    glGenVertexArrays(1, &quad_VertexArrayID);
    glBindVertexArray(quad_VertexArrayID);

    static const GLfloat g_quad_vertex_buffer_data[] = {
        -1.0f, -1.0f, 0.0f, 1.0f, -1.0f, 0.0f, -1.0f, 1.0f, 0.0f, -1.0f, 1.0f, 0.0f, 1.0f, -1.0f, 0.0f, 1.0f, 1.0f, 0.0f,
    }; GLuint vao, textureID;
    // Create one OpenGL texture
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glGenTextures(1, &textureID);

    GLuint quad_vertexbuffer;
    glGenBuffers(1, &quad_vertexbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, quad_vertexbuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_quad_vertex_buffer_data), g_quad_vertex_buffer_data, GL_STATIC_DRAW);

    // Create and compile our GLSL program from the shaders
    GLint program = createProgram(vtxShader, frgShader);

    glUseProgram(program);
    GLint texID = glGetUniformLocation(program, "renderedTexture");
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, textureID);
    glUniform1i(texID, 0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, quad_vertexbuffer);
    glVertexAttribPointer(0,        // attribute 0. No particular reason for 0, but
                                    // must match the layout in the shader.
        3,        // size
        GL_FLOAT, // type
        GL_FALSE, // normalized?
        0,        // stride
        (void*)0 // array buffer offset
    );
    glUseProgram(program);

    glUseProgram(0);
    glBindVertexArray(0);

    static int32_t* data;

    static int32_t oldWidth = -1, oldHeight = -1;
    if (oldWidth != screenWidth || oldHeight != screenHeight) {
        oldWidth = screenWidth;
        oldHeight = screenHeight;
        if (data != nullptr) {
            free(data);
        }
        data = new int32_t[screenWidth * screenHeight];
        for (int32_t j = 0; j < screenHeight; ++j) {
            for (int32_t i = 0; i < screenWidth; ++i) {
                data[j * screenWidth + i] = rand();
            }
            //break;
        }
        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, screenWidth, screenHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    }
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);

        if (oldWidth != screenWidth || oldHeight != screenHeight) {
            oldWidth = screenWidth;
            oldHeight = screenHeight;
            if (data != nullptr) {
                free(data);
            }
            data = new int32_t[screenWidth * screenHeight];
            for (int32_t j = 0; j < screenHeight; ++j) {
                for (int32_t i = 0; i < screenWidth; ++i) {
                    data[j * screenWidth + i] = rand();
                }
                //break;
            }
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, textureID);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, screenWidth, screenHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        }


        auto dx = xmax - xmin;
        auto dy = ymax - ymin;
        char* dataC = (char*)data;
        for (int32_t j = 0; j < screenHeight; ++j) {
            for (int32_t i = 0; i < screenWidth; ++i) {
                auto xrel = ((double)i) / ((double)screenWidth);
                auto yrel = ((double)j) / ((double)screenHeight);

                auto xpos = xmin + xrel * dx;
                auto ypos = ymin + yrel * dy;

                auto ddx = data_xmax - data_xmin;
                auto ddy = data_ymax - data_ymin;
                if (xpos >= data_xmin && xpos <= data_xmax && ypos >= data_ymin && ypos <= data_ymax) {
					auto relx = std::clamp((xpos - data_xmin) / (data_xmax - data_xmin), 0., 1.);
                    auto rely = std::clamp((ypos - data_ymin) / (data_ymax - data_ymin), 0., 1.);
                    auto l = (int32_t) ::floor(relx * ((double)xLim - 1.));
                    auto r = (int32_t) ::ceil(relx * ((double)xLim - 1.));
                    auto b = (int32_t) ::floor(rely * ((double)yLim - 1.));
                    auto t = (int32_t) ::ceil(rely * ((double)yLim - 1.));

                    auto alpha = relx * ((double)xLim - 1.) - ::floor(relx * ((double)xLim - 1.));
                    auto beta = rely * ((double)yLim - 1.) - ::floor(rely * ((double)yLim - 1.));

                    auto int_to_v4 = [](int32_t v) {
                        char* raw = reinterpret_cast<char*>(&v);
                        Eigen::Vector4d vec(raw[0], raw[1], raw[2], raw[3]);
                        return vec;
                    };
                    auto v4_to_int = [](Eigen::Vector4d vec) {
                        int32_t v = 0;
                        char* raw = reinterpret_cast<char*>(&v);
                        raw[0] = vec.x();
                        raw[1] = vec.y();
                        raw[2] = vec.z();
                        raw[3] = vec.w();
                        return v;
                    };

                    auto vbl = int_to_v4(imageBuffer[b * xLim + l]);
                    //auto vbr = int_to_v4(imageBuffer[b * xLim + r]);
                    //auto vtl = int_to_v4(imageBuffer[t * xLim + l]);
                    //auto vtr = int_to_v4(imageBuffer[t * xLim + r]);

                    //auto vb = vbr * alpha + vbl * (1. - alpha);
                    //auto vt = vtr * alpha + vtl * (1. - alpha);

                    //auto v = vt * beta + vb * (1. - beta);
                    data[j * screenWidth + i] = v4_to_int(vbl);
                    //data[j * screenWidth + i] = imageBuffer[b * xLim + l];

                    //dataC[(j * screenWidth + i) * 4 + 0] = xrel * 255.f;
                    //dataC[(j * screenWidth + i) * 4 + 1] = yrel * 255.f;
                    //dataC[(j * screenWidth + i) * 4 + 2] = 0;
                    dataC[(j * screenWidth + i) * 4 + 3] = 255;
                }
                else {
                    dataC[(j * screenWidth + i) * 4 + 0] = 0;
                    dataC[(j * screenWidth + i) * 4 + 1] = 0;
                    dataC[(j * screenWidth + i) * 4 + 2] = 0;
                    dataC[(j * screenWidth + i) * 4 + 3] = 255;
                }
            }
        }
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, screenWidth, screenHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);



        glUseProgram(program);
        glBindVertexArray(vao);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, textureID);
        //glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, screenWidth, screenHeight, GL_RGBA, GL_UNSIGNED_BYTE, data);

        glDrawArrays(GL_TRIANGLES, 0, 6); // 2*3 indices starting at 0 -> 2 triangles
        glBindVertexArray(0);
        glUseProgram(0);

        ImGui::NewFrame();
        bool show_parameter_window = true;
        //ParameterManager::instance().buildImguiWindow(&show_parameter_window);
        OSD();
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable) {
            GLFWwindow* backup_current_context = glfwGetCurrentContext();
            ImGui::UpdatePlatformWindows();
            ImGui::RenderPlatformWindowsDefault();
            glfwMakeContextCurrent(backup_current_context);
        }
        glfwSwapBuffers(window);
    }

    return 0;
}