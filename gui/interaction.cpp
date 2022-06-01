#include "glui.h"
#include <imgui/imgui.h>
#include <iostream>

void GUI::keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	static bool emit = false;
	if (action != GLFW_PRESS) { 
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
	}
}
bool trackingLeft = false;
bool trackingRight = false;

void GUI::mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
	auto [minDomain, maxDomain] = getDomain();


	// if (action == GLFW_PRESS) {
	// 	if (button == GLFW_MOUSE_BUTTON_LEFT) {
	// 		trackingLeft = true;
	// 		double x, y; glfwGetCursorPos(window, &x, &y);
	// 		orig = vec(x / (scalar)screenWidth * domainWidth, (screenHeight - y) / (scalar)screenHeight * domainHeight);
	// 	}
	// 	if (button == GLFW_MOUSE_BUTTON_RIGHT) {
	// 		trackingRight = true;
	// 		double x, y; glfwGetCursorPos(window, &x, &y);
	// 		target = vec(x / (scalar)screenWidth * domainWidth, (screenHeight - y) / (scalar)screenHeight * domainHeight);
	// 	}
	// }
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
	// if (trackingLeft) {
	// 	orig = vec(xpos / (scalar)screenWidth * domainWidth, (screenHeight - ypos) / (scalar)screenHeight * domainHeight);
	// }
	// if (trackingRight) {
	// 	//std::cout << xpos << " x " << ypos << " -> " << xpos / (scalar)screenWidth * domainWidth << " x " << ypos / (scalar)screenHeight * domainHeight << std::endl;
	// 	target = vec(xpos / (scalar)screenWidth * domainWidth, (screenHeight - ypos) / (scalar)screenHeight * domainHeight);
	// }
	m_cursorPosition = vec(xpos, ypos);


}
void GUI::scrollCallback(GLFWwindow* window, double xoffset, double yoffset) { }
