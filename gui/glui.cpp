#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"
#include "glui.h"

GUI& GUI::instance() {
	static GUI instance;
	return instance;
}
GUI::GUI() {
    //ParameterManager::instance().newParameter("camera.perspective", glm::mat4());
    //ParameterManager::instance().newParameter("camera.view", glm::mat4());
    //ParameterManager::instance().newParameter("simTime", float{ 0.f }, { .constant = true });
    //ParameterManager::instance().newParameter("renderTime", float{ 0.f }, { .constant = true });
}
void GUI::uiFunctions() {
    if (show_parameter_window) ParameterManager::instance().buildImguiWindow(&show_parameter_window);
    if (show_timer_window) TimerWindow(&show_timer_window);
    static bool show_ray_window = true;
    RayWindow(&show_ray_window);
    OSD();
    ImGui::Render();
}
void GUI::quit() {
    shouldStop = true;
    if (m_ffmpegPipe != nullptr)
#ifdef WIN32
        _pclose(m_ffmpegPipe);
#else
        pclose(m_ffmpegPipe);
#endif
    if (summaryFileOpen)
        summaryFile.close();
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    render_lock.unlock();
}
