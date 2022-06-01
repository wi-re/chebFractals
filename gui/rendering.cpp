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
        // auto dt = 1.0 / 400.0;
        // static auto lastTime = ParameterManager::instance().get<double>("sim.time") - dt;
        // auto currTime = ParameterManager::instance().get<double>("sim.time");
        // if (currTime > lastTime + dt && m_ffmpegPipe != nullptr) {
        //     lastTime = currTime;
        //     static int32_t* buffer = new int32_t[screenWidth * screenHeight];
        //     glReadPixels(0, 0, screenWidth, screenHeight, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
        //     fwrite(buffer, sizeof(int) * screenWidth * screenHeight, 1, m_ffmpegPipe);

        // }

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
        pclose(m_ffmpegPipe);
    quit();
}
