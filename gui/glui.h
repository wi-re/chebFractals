#pragma once
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "imgui/imgui.h"
#include <fsVis/utility.h>
#include <fsVis/optimizers.h>
#include <mutex>
#include <vector>
#include <tools/ParameterManager.h>
#include <glm/glm.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

template <typename T = float>
auto bind(GLuint program, GLuint location, T value)
{
    using Ty = std::decay_t<T>;
    glUseProgram(program);
    if constexpr (std::is_same_v<Ty, int32_t>)
        glUniform1i(location, value);
    else if constexpr (std::is_same_v<Ty, bool>)
        glUniform1i(location, value ? 1 : 0);
    else if constexpr (std::is_same_v<Ty, uint32_t>)
        glUniform1ui(location, value);
    else if constexpr (std::is_same_v<Ty, float>)
        glUniform1f(location, value);
    else if constexpr (std::is_same_v<Ty, glm::mat4>)
        glUniformMatrix4fv(location, 1, false, &value[0][0]);
    glUseProgram(0);
}

struct gl_uniform_base
{
    std::vector<std::pair<GLuint, GLuint>> programs;

    virtual void update() = 0;
    virtual void add_uniform(GLuint program) = 0;
    virtual void add_uniform(GLuint program, GLuint location) = 0;
};

class GUI
{
    GUI();
    bool show_parameter_window = true;
    GLFWwindow *window;
    ImVec4 clear_color;
    bool exportFlag = false;
    bool m_showText = true;
    bool m_showInfo = true;
    vec m_cursorPosition;
    bool pickerActive;
    std::map<std::string, gl_uniform_base *> m_parameterMappings;
    std::map<std::string, gl_uniform_base *> m_uniformMappings;
    std::map<std::string, std::vector<std::variant<int32_t, float, std::vector<int32_t>, std::vector<float>>>> m_renderData;

    bool prettyRender = false;
    bool writeFlag = false;
    FILE *m_ffmpegPipe = nullptr;

public:
    bool simulationRunning = false;
    bool shouldStop = false;
    std::mutex render_lock;
    static GUI &instance();

    void renderFunctions();
    void uiFunctions();

    void RayWindow(bool *p_open);

    void initGL();
    void renderLoop();
    void quit();
    void initParameters(int argc, char *argv[]);
    void initGVDB();
    void initSimulation();

    void OSD();

    void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods);
    void cursorPositionCallback(GLFWwindow *window, double xpos, double ypos);
    void mouseButtonCallback(GLFWwindow *window, int button, int action, int mods);
    void scrollCallback(GLFWwindow *window, double xoffset, double yoffset);
};
