#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"
#include <glad/glad.h> 
#include "glui.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <tools/stb_image_write.h>
#define STB_IMAGE_READ_IMPLEMENTATION
#include <tools/stb_image.h>
#include <simulation/SPH.h>

GLuint shader_programme;
GLuint vao = 0;
void GUI::renderFunctions() {
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glEnable(GL_MULTISAMPLE);
    glClearColor(0.2f, 0.2f, 0.2f, 1.f);
    glClearColor(1.0f, 1.0f, 1.0f, 1.f);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHT0);
    glDisable(GL_LIGHTING);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glDisable(GL_COLOR_MATERIAL);
    glFlush();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

        render();

    static int32_t i = 0;
    // ACTUAL RENDERING CALL

    glBindVertexArray(0);
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glPopAttrib();
}
void GUI::renderLoop() {
    static auto& xmin = ParameterManager::instance().get<scalar>("domain.xmin");
    static auto& xmax = ParameterManager::instance().get<scalar>("domain.xmax");
    static auto& ymin = ParameterManager::instance().get<scalar>("domain.ymin");
    static auto& ymax = ParameterManager::instance().get<scalar>("domain.ymax");
    //xmin = 0x1.e83a227ae4c3ep-3; ymin = -0x1.9b64feac755ffp-26; xmax = 0x1.e83a227ae4ddap-3;ymax =  -0x1.9b64f1b139c55p-26;
    

    show_demo_window = true;
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    while (!shouldStop && !glfwWindowShouldClose(window)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);

        renderFunctions();
        {
            if(exportFlag){
                exportFlag = false;
		
        const std::chrono::time_point<std::chrono::system_clock> now =
            std::chrono::system_clock::now();
        const std::time_t t_c = std::chrono::system_clock::to_time_t(now);
        auto time = std::time(nullptr);
        std::stringstream ss;
        char buf[256];
        ss << std::put_time(std::localtime(&t_c), "%F_%T") 
		<< "_" << ParameterManager::instance().get<std::string>("field.method")
		<< ".png"; // ISO 8601 without timezone information.
        auto s = ss.str();
        std::replace(s.begin(), s.end(), ':', '-');
		
            static int32_t* buffer = new int32_t[screenWidth * screenHeight];
            char* buffC = (char*) buffer;
            glReadPixels(0, 0, screenWidth, screenHeight, GL_RGBA, GL_UNSIGNED_BYTE, buffer);

            static char* bufferC = new char[screenWidth * screenHeight * 3];
            for(int32_t i = 0; i < screenWidth; i++){
                for(int32_t j = 0; j < screenHeight; j++){
                    bufferC[(i + screenWidth * ( screenHeight -j - 1))*3 + 0] = buffC[(i + screenWidth * j ) * 4 + 0];
                    bufferC[(i + screenWidth * ( screenHeight -j - 1))*3 + 1] = buffC[(i + screenWidth * j ) * 4 + 1];
                    bufferC[(i + screenWidth * ( screenHeight -j - 1))*3 + 2] = buffC[(i + screenWidth * j ) * 4 + 2];
                }
            }


if (!stbi_write_png(s.c_str(), screenWidth, screenHeight,3, bufferC, screenWidth*3)) {
  std::cerr << "ERROR: could not write image to " << s << std::endl;
}
            }
        }

#define WATCH(type, ns, id)\
  static type id = ParameterManager::instance().get<type>(std::string(#ns) +"." + #id);\
  if (ParameterManager::instance().get<type>(std::string(#ns) +"." + #id) != id) {\
      std::cout << "Parameter " << std::string(#ns) << "." << std::string(#id) << " changed from " << id << " to " << ParameterManager::instance().get<type>(std::string(#ns) +"." + #id) << std::endl;\
      id = ParameterManager::instance().get<type>(std::string(#ns) +"." + #id);\
      dirty = true;\
  }
        bool dirty = false;
        WATCH(scalar, domain, xmin);
        WATCH(scalar, domain, xmax);
        WATCH(scalar, domain, ymin);
        WATCH(scalar, domain, ymax);


         if (m_ffmpegPipe != nullptr && dirty) {
             static int32_t* buffer = new int32_t[screenWidth * screenHeight];
             glReadPixels(0, 0, screenWidth, screenHeight, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
             fwrite(buffer, sizeof(int) * screenWidth * screenHeight, 1, m_ffmpegPipe);

         }

        ImGui::NewFrame();
        //ImGui::ShowDemoWindow(&show_demo_window);
        uiFunctions();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable)
        {
            GLFWwindow* backup_current_context = glfwGetCurrentContext();
            ImGui::UpdatePlatformWindows();
            ImGui::RenderPlatformWindowsDefault();
            glfwMakeContextCurrent(backup_current_context);
        }
        glfwSwapBuffers(window);

    }
    if (m_ffmpegPipe != nullptr)
#ifdef WIN32
        _pclose(m_ffmpegPipe);
#else
        pclose(m_ffmpegPipe);
#endif
    quit();
}
