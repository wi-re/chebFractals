#include "glui.h"
void GUI::OSD(){
    if (m_showText) {
        const float DISTANCE = 10.0f;
        static int corner = 0;
        ImGuiIO& io = ImGui::GetIO();
        if (corner != -1)
        {
            ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImVec2 work_area_pos = viewport->WorkPos;   // Instead of using viewport->Pos we use GetWorkPos() to avoid menu bars, if any!
            ImVec2 work_area_size = viewport->WorkSize;
            ImVec2 window_pos = ImVec2((corner & 1) ? (work_area_pos.x + work_area_size.x - DISTANCE) : (work_area_pos.x + DISTANCE), (corner & 2) ? (work_area_pos.y + work_area_size.y - DISTANCE) : (work_area_pos.y + DISTANCE));
            ImVec2 window_pos_pivot = ImVec2((corner & 1) ? 1.0f : 0.0f, (corner & 2) ? 1.0f : 0.0f);
            ImGui::SetNextWindowPos(window_pos, ImGuiCond_Always, window_pos_pivot);
            ImGui::SetNextWindowViewport(viewport->ID);
        }
        ImGui::SetNextWindowBgAlpha(0.35f); // Transparent background
        if (ImGui::Begin("Stats overlay", &m_showText, (corner != -1 ? ImGuiWindowFlags_NoMove : 0) | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav))
        {
            auto addParameter = [&](auto paramName) {
                auto& param = ParameterManager::instance().getParameter(paramName);
                ImGui::PushID(param.identifier.c_str());
                ParameterManager::instance().uiFunctions[param.type](param);
                ImGui::PopID();
            };

            //addParameter("color_map.buffer");
            addParameter("field.min");
            addParameter("field.max");

        }
        ImGui::End();
    }
    if(m_showInfo){
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
        //std::cout << "[ " << minDomain.x() << " : " << minDomain.y() << " ] -> [ " << maxDomain.x() << " : " << maxDomain.y() << " ]\n"; 
        //std::cout << m_cursorPosition.x() << " : " << m_cursorPosition.y() << " -> " << x << " : " << y << std::endl;

        vec position(x, y);

        const float DISTANCE = 10.0f;
        static int corner = 3;
        ImGuiIO& io = ImGui::GetIO();
        if (corner != -1)
        {
            ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImVec2 work_area_pos = viewport->WorkPos;   // Instead of using viewport->Pos we use GetWorkPos() to avoid menu bars, if any!
            ImVec2 work_area_size = viewport->WorkSize;
            ImVec2 window_pos = ImVec2((corner & 1) ? (work_area_pos.x + work_area_size.x - DISTANCE) : (work_area_pos.x + DISTANCE), (corner & 2) ? (work_area_pos.y + work_area_size.y - DISTANCE) : (work_area_pos.y + DISTANCE));
            ImVec2 window_pos_pivot = ImVec2((corner & 1) ? 1.0f : 0.0f, (corner & 2) ? 1.0f : 0.0f);
            ImGui::SetNextWindowPos(window_pos, ImGuiCond_Always, window_pos_pivot);
            ImGui::SetNextWindowViewport(viewport->ID);
        }
        ImGui::SetNextWindowBgAlpha(0.35f); // Transparent background
        if (ImGui::Begin("Particle Info", &m_showInfo, (corner != -1 ? ImGuiWindowFlags_NoMove : 0) | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav))
        {
            auto addScalar = [&](auto name, scalar value, auto description) {
                float vcp = value;
                auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
                ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
                ImGui::DragFloat(name, &vcp, 0, vcp, vcp,"%.6f");
                ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
                if(ImGui::IsItemHovered())
                    ImGui::SetTooltip("%s",description);
                return;
            };
            struct float2{float x,y;};
            auto addScalar2 = [&](auto name, vec value, auto description) {
                float2 vcp(value.x(), value.y());
                auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
                ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
                ImGui::DragFloat2(name, &vcp.x, 0, 0, 0, "%.6f");
                ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("%s",description);
                return;
            };
            auto addInteger = [&](auto name, int32_t value, auto description) {
                int32_t vcp = value;
                auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
                ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
                ImGui::DragInt(name, &vcp, 0, vcp, vcp);
                ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("%s",description);
                return;
            };



            addScalar2("Cursor        ", position, "");
            {                
    auto ny = ParameterManager::instance().get<int32_t>("field.ny");
    auto nx = ParameterManager::instance().get<int32_t>("field.nx");
                auto x = std::clamp(m_cursorPosition.x(), 0.0, (scalar)screenWidth) / (scalar)screenWidth;
                auto y = 1. - std::clamp(m_cursorPosition.y(), 0.0, (scalar)screenHeight) / (scalar)screenHeight;

                int32_t xi = std::clamp((int32_t)::floor(x * dataWidth),0,nx-1);
                int32_t yi = std::clamp((int32_t)::floor(y * dataHeight),0,ny-1);

                scalar value = scalarFieldData[xi + yi * dataWidth];
                addScalar("Scalar Field", value, "");
            }


        }
        ImGui::End();
    }
}
