#pragma once
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "imgui/imgui.h"
#include <simulation/SPH.h>
#include <mutex>
#include <vector>
#include <tools/ParameterManager.h>
#include <glm/glm.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

template<typename T = float>
auto bind(GLuint program, GLuint location, T value) {
    using Ty = std::decay_t<T>;
    glUseProgram(program);
    if constexpr (std::is_same_v<Ty, int32_t>)
        glUniform1i(location, value);
    else if constexpr (std::is_same_v<Ty, bool>)
        glUniform1i(location, value ? 1: 0);
    else if constexpr (std::is_same_v<Ty, uint32_t>)
        glUniform1ui(location, value);
    else if constexpr (std::is_same_v<Ty, float>)
        glUniform1f(location, value);
    // else if constexpr (std::is_same_v<Ty, float2>)
    //     glUniform2f(location, value.x, value.y);
    // else if constexpr (std::is_same_v<Ty, float3>)
    //     glUniform3f(location, value.x, value.y, value.z);
    // else if constexpr (std::is_same_v<Ty, float4>)
    //     glUniform4f(location, value.x, value.y, value.z, value.w);
    // else if constexpr (std::is_same_v<Ty, int3>)
    //     glUniform3i(location, value.x, value.y, value.z);
    // else if constexpr (std::is_same_v<Ty, int4>)
    //     glUniform4i(location, value.x, value.y, value.z, value.w);
    else if constexpr (std::is_same_v<Ty, glm::mat4>)
        glUniformMatrix4fv(location, 1, false, &value[0][0]);
    // else
        // std::enable_if_t<false>();
    glUseProgram(0);
}

struct gl_uniform_base {
    std::vector<std::pair<GLuint, GLuint>> programs;

    virtual void update() = 0;
    virtual void add_uniform(GLuint program) = 0;
    virtual void add_uniform(GLuint program, GLuint location) = 0;
};
// template <typename T> struct gl_uniform_custom : public gl_uniform_base {
//     T* ptr;
//     std::function<bool()> update_fn;
//     std::string variableName;

//     gl_uniform_custom(T* p, std::string identifier)
//         : ptr(p), variableName(identifier) {
//         update_fn = [this]() {
//             static T old = *ptr;
//             static auto old_size = programs.size();

//             if (old != *ptr || old_size != programs.size()) {
//                 old = *ptr;
//                 old_size = programs.size();
//                 return true;
//             }
//             return false;
//         };
//         update_fn();
//     }
//     virtual void update() {
//         for (auto p : programs) {
//             auto [program, uniform] = p;
//             bind(program, uniform, *ptr);
//             //program->bind();
//             //std::decay_t<decltype(*ptr)> val = *ptr;
//             //program->setUniformValue(uniform, convertToQt(val));
//             //program->release();
//         }
//     }
//     virtual void add_uniform(GLuint prog, GLuint location) {
//         programs.emplace_back(prog, location);
//     }
//     virtual void add_uniform(GLuint prog) {
//         auto identifier = glGetUniformLocation(prog, variableName.c_str());
//         //auto identifier = prog->uniformLocation(variableName.c_str());
//         if (identifier == -1)
//             return;
//         //std::cout << "Found uniform for program " << prog << ": " << variableName.c_str() << std::endl;
//         programs.emplace_back(prog, identifier);
//     }
// };
// template<typename Ty, Ty ident> struct gl_uniform_parameter : public gl_uniform_base {
//     using T = std::decay_t<decltype(get<ident>())>;
//     std::string glIdentifier;
//     gl_uniform_parameter() {
//         auto nsid = getIdentifier<ident>();
//         auto ns = nsid.first;
//         auto id = nsid.second;
//         if (ParameterManager::instance().isAmbiguous(id) == true) {
//             glIdentifier = ns + "_" + id;
//         }
//         else
//             glIdentifier = id;
//     }
//     virtual void update() {
//         for (auto p : programs) {
//             auto [program, uniform] = p;
//             glUseProgram(program);
//             //program->bind();
//             bind(program, uniform, get<ident>());
//             //std::decay_t<decltype(*Ty::ptr)> val = *Ty::ptr;
//             //glUniform
//             //    program->setUniformValue(uniform, convertToQt(val));
//             program->release();
//         }
//     }
//     virtual void add_uniform(GLuint prog, GLuint location) {
//         programs.emplace_back(prog, location);
//     }
//     virtual void add_uniform(GLuint prog) {
//         auto identifier = glGetUniformLocation(prog, glIdentifier.c_str());
//         //auto identifier = prog->uniformLocation(variableName.c_str());
//         if (identifier == -1)
//             return;
//         //std::cout << "Found uniform for program " << prog << ": " << variableName.c_str() << std::endl;
//         programs.emplace_back(prog, identifier);
//     }
// };
// #if __has_include(<cuda_runtime.h>)
// struct cuda_buffer_base {
//     cudaGraphicsResource* resource = nullptr;
//     GLuint VBO = UINT_MAX;

//     virtual ~cuda_buffer_base() {
//         if (resource != nullptr) {
//             cudaGraphicsUnregisterResource(resource);
//             resource = nullptr;
//         }
//         if (VBO != UINT_MAX) {
//             glDeleteBuffers(1, &VBO);
//             VBO = UINT_MAX;
//         }
//     }
//     virtual void update() = 0;
//     virtual GLint bindAttribute(GLuint attribute) = 0;
//     virtual GLint bindProgram(GLuint prog) = 0;
// };
// template <typename info>
// struct cuda_buffer : public cuda_buffer_base {
//     using T = typename info::type;
//     bool initialized = false;
//     T* dptr;
//     size_t num_bytes;

//     std::function<bool()> update_fn;

//     void initMemory() {
//         //initialized = false;
//         //return;
//         uint32_t max_numptcls = get<parameters::simulation_settings::maxNumptcls>();
//         glGenBuffers(1, &VBO);
//         glBindBuffer(GL_ARRAY_BUFFER, VBO);
//         glBufferData(GL_ARRAY_BUFFER, max_numptcls * sizeof(T), 0, GL_STREAM_DRAW);
//         cudaGraphicsGLRegisterBuffer(&resource, VBO,
//             cudaGraphicsMapFlagsWriteDiscard);
//         glBindBuffer(GL_ARRAY_BUFFER, 0);
//         initialized = true;
//     }
//     virtual ~cuda_buffer() {
//         if (resource != nullptr) {
//             cudaGraphicsUnregisterResource(resource);
//             resource = nullptr;
//         }
//         if (VBO != UINT_MAX) {
//             glDeleteBuffers(1, &VBO);
//             VBO = UINT_MAX;
//         }
//     }
//     cuda_buffer(std::function<bool()> valid = []() {
//         static auto old_frame = INT_MAX;
//         if (old_frame != get<parameters::internal::frame>()) {
//             old_frame = get<parameters::internal::frame>();
//             return true;
//         }
//         return false;
//         }) : update_fn(valid) {
//         //initializeOpenGLFunctions();
//     }
//         virtual void update() {
//             if (!update_fn())
//                 return;
//             if (!initialized)
//                 return;
//             cudaGraphicsMapResources(1, &resource, 0);
//             cudaGraphicsResourceGetMappedPointer((void**)&dptr, &num_bytes, resource);
//             if (get<parameters::internal::target>() != launch_config::device)
//                 cudaMemcpy(dptr, info::ptr, info::alloc_size, cudaMemcpyHostToDevice);
//             else
//                 cudaMemcpy(dptr, info::ptr, info::alloc_size, cudaMemcpyDeviceToDevice);
//             cudaGraphicsUnmapResources(1, &resource, 0);
//         }
//         virtual GLint bindAttribute(GLuint attribute) {
//             if (!initialized)
//                 initMemory();
//             //return -1;
//             glBindBuffer(GL_ARRAY_BUFFER, VBO);
//             glEnableVertexAttribArray(attribute);
//             glVertexAttribPointer(attribute, sizeof(T) / sizeof(float), GL_FLOAT,
//                 GL_FALSE, 0, NULL);
//             glVertexAttribDivisor(attribute, 1);
//             glBindBuffer(GL_ARRAY_BUFFER, 0);
//             return attribute;
//         }
//         virtual GLint bindProgram(GLuint prog) {
//            // glUseProgram(prog);
//             std::string name = info::qualifiedName;
//             std::vector<std::string> SplitVec;
//             boost::algorithm::split(SplitVec, name, boost::algorithm::is_any_of("."), boost::algorithm::token_compress_on);
//             auto attribute = glGetAttribLocation(prog, SplitVec[1].c_str());
//             //auto attribute = prog->attributeLocation(SplitVec[1].c_str());
//             if (attribute == -1) {
//                 attribute = glGetAttribLocation(prog, name.c_str());
//                 //attribute = prog->attributeLocation(name.c_str());
//                 if (attribute == -1) {
//                    // glUseProgram(0);
//                     return-1;
//                 }
//                 //programs.emplace_back(prog, attribute);
//                 //return;
//             }
//             //auto attribute = prog->attributeLocation(info::variableName);
//             if (attribute == -1) {
//                 return -1;
//               //  glUseProgram(0);
//             }
//             if (!initialized)
//                 initMemory();
//             //std::cout << "Found attribute " << name << " for program " << prog << " @ " << attribute << std::endl;
//             //return -1;
//             glBindBuffer(GL_ARRAY_BUFFER, VBO);
//             glEnableVertexAttribArray(attribute);
//             glVertexAttribPointer(attribute, sizeof(T) / sizeof(float), GL_FLOAT,
//                 GL_FALSE, 0, NULL);
//             glVertexAttribDivisor(attribute, 1);
//             glBindBuffer(GL_ARRAY_BUFFER, 0);
//            // glUseProgram(0);
//             return attribute;
//         }
// };
// #endif
class GUI {
	GUI();
	bool show_demo_window = true;
	bool show_parameter_window = true;
	bool show_timer_window = true;
	bool show_log_window = true;
	GLFWwindow* window;
	ImVec4 clear_color;

	bool m_showText = true;
    bool m_showInfo = true;
    vec m_cursorPosition;
	bool pickerActive;
    	std::map<std::string, gl_uniform_base *> m_parameterMappings;
    	std::map<std::string, gl_uniform_base *> m_uniformMappings;
        std::map<std::string, std::vector<std::variant<int32_t, float, std::vector<int32_t>, std::vector<float>>>> m_renderData;

        bool prettyRender = false;
        bool writeFlag = false;
        FILE* m_ffmpegPipe = nullptr;
public:
    bool simulationRunning = false;
	bool shouldStop = false;
	std::mutex render_lock;
	static  GUI& instance();

	void renderFunctions();
	void uiFunctions();

	void TimerWindow(bool* p_open);
    void RayWindow(bool* p_open);

	void initGL();
	void renderLoop();
	void quit();
	void initParameters(int argc, char* argv[]);
	void initGVDB();
	void initSimulation();

    void OSD();

	void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
	void cursorPositionCallback(GLFWwindow* window, double xpos, double ypos);
	void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
	void scrollCallback(GLFWwindow* window, double xoffset, double yoffset);
};

