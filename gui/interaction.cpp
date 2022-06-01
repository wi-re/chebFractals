#include "glui.h"
#include <imgui/imgui.h>
#include <iostream>

void GUI::keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
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
		case GLFW_KEY_LEFT: xmin -= dx * (scaling - 1.); xmax -= dx * (scaling - 1.); break;
		case GLFW_KEY_RIGHT: xmin += dx * (scaling - 1.); xmax += dx * (scaling - 1.); break;
		case GLFW_KEY_UP: ymin += dy * (scaling - 1.); ymax += dy * (scaling - 1.); break;
		case GLFW_KEY_DOWN: ymin -= dy * (scaling - 1.); ymax -= dy * (scaling - 1.); break;
		case GLFW_KEY_APOSTROPHE: xmin = xcenter - dx / 2. * scaling; xmax = xcenter + dx / 2. * scaling; ymin = ycenter - dy / 2. * scaling; ymax = ycenter + dy / 2. * scaling; break;
		case GLFW_KEY_RIGHT_BRACKET: xmin = xcenter - dx / 2. / scaling; xmax = xcenter + dx / 2. / scaling; ymin = ycenter - dy / 2. / scaling; ymax = ycenter + dy / 2. / scaling; break;
		}

		if (key == GLFW_KEY_E && action == GLFW_RELEASE)
			emit = false;
		
		return; }
	switch (key) {
	case GLFW_KEY_H:
		m_showText = !m_showText;
		break;
	case GLFW_KEY_G:
	{
		auto& grid = ParameterManager::instance().get<bool>("render.showGrid");
		grid = !grid;
		break;
	}
	case GLFW_KEY_1: ParameterManager::instance().get<int32_t>("colorMap.map") = 0; break;
	case GLFW_KEY_2: ParameterManager::instance().get<int32_t>("colorMap.map") = 1; break;
	case GLFW_KEY_3: ParameterManager::instance().get<int32_t>("colorMap.map") = 2; break;
	case GLFW_KEY_T: ParameterManager::instance().get<bool>("colorMap.auto") = !ParameterManager::instance().get<bool>("colorMap.auto"); break;
	case GLFW_KEY_F: ParameterManager::instance().get<bool>("field.render") = !ParameterManager::instance().get<bool>("field.render"); break;
	case GLFW_KEY_M: ParameterManager::instance().get<bool>("marching.render") = !ParameterManager::instance().get<bool>("marching.render"); break;
	case GLFW_KEY_N: ParameterManager::instance().get<bool>("marching.solid") = !ParameterManager::instance().get<bool>("marching.solid"); break;
	case GLFW_KEY_O: ParameterManager::instance().get<bool>("ptcl.render") = !ParameterManager::instance().get<bool>("ptcl.render"); break;
	case GLFW_KEY_E: emit=true; break;
	case GLFW_KEY_L: std::cout << std::setprecision(16) << std::hexfloat << xmin << " " << ymin << " " << xmax << " " << ymax << std::endl << std::defaultfloat; break;


	}

	{


	}

}
bool trackingLeft = false;
bool trackingRight = false;

void GUI::mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
auto [minDomain, maxDomain] = getDomain();


        auto x = std::clamp(m_cursorPosition.x(), 0.0, (scalar)screenWidth) / (scalar)screenWidth;
        auto y = std::clamp(m_cursorPosition.y(), 0.0, (scalar)screenHeight) / (scalar)screenHeight;
        if (y > 1.0)
            y -= 1.0;
        y = 1.0 - y;
        x *= maxDomain.x() - minDomain.x();
        y *= maxDomain.y() - minDomain.y();
        x+=minDomain.x();
        y+=minDomain.y();

	if (action == GLFW_PRESS) {
		if (button == GLFW_MOUSE_BUTTON_LEFT) {
			trackingLeft = true;
			ParameterManager::instance().get<scalar>("path.x") = x;
			ParameterManager::instance().get<scalar>("path.y") = y;
		}
	// 	if (button == GLFW_MOUSE_BUTTON_RIGHT) {
	// 		trackingRight = true;
	// 		double x, y; glfwGetCursorPos(window, &x, &y);
	// 		target = vec(x / (scalar)screenWidth * domainWidth, (screenHeight - y) / (scalar)screenHeight * domainHeight);
	// 	}
	}
	if(action == GLFW_RELEASE) {
		if (button == GLFW_MOUSE_BUTTON_LEFT)
			trackingLeft = false;
		if (button == GLFW_MOUSE_BUTTON_RIGHT)
			trackingRight = false;
	}
}

//#include <glrender/glparticleIndexRender/particleIndexRender.h>
#include <imgui/imgui_internal.h>
void GUI::cursorPositionCallback(GLFWwindow* window, double xpos, double ypos) {
	auto GImGui = ImGui::GetCurrentContext();
	// if (trackingRight) {
	// 	//std::cout << xpos << " x " << ypos << " -> " << xpos / (scalar)screenWidth * domainWidth << " x " << ypos / (scalar)screenHeight * domainHeight << std::endl;
	// 	target = vec(xpos / (scalar)screenWidth * domainWidth, (screenHeight - ypos) / (scalar)screenHeight * domainHeight);
	// }
	m_cursorPosition = vec(xpos, ypos);

auto [minDomain, maxDomain] = getDomain();


        auto x = std::clamp(m_cursorPosition.x(), 0.0, (scalar)screenWidth) / (scalar)screenWidth;
        auto y = std::clamp(m_cursorPosition.y(), 0.0, (scalar)screenHeight) / (scalar)screenHeight;
        if (y > 1.0)
            y -= 1.0;
        y = 1.0 - y;
        x *= maxDomain.x() - minDomain.x();
        y *= maxDomain.y() - minDomain.y();
        x+=minDomain.x();
        y+=minDomain.y();
	if (trackingLeft) {
		
			ParameterManager::instance().get<scalar>("path.x") = x;
			ParameterManager::instance().get<scalar>("path.y") = y;
	}



}
void GUI::scrollCallback(GLFWwindow* window, double xoffset, double yoffset) { }
