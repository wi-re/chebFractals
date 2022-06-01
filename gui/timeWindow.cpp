#include <glad/glad.h> 
#include "glui.h"
#include <tools/timer.h>

ImVec4 hexToCol(Color hexc) {
    int64_t hex = (int64_t)hexc;
    auto r = (hex & 0x00FF0000ul) >> 16;
    auto g = (hex & 0x0000FF00ul) >> 8;
    auto b = (hex & 0x000000FFul) >> 0;
    auto a = (hex & 0xFF000000ul) >> 24;
    ImVec4 col(
        (float)(r) / 255.f,
        (float)(g) / 255.f,
        (float)(b) / 255.f,
        (float)(a) / 255.f);
    return col;
}

void GUI::TimerWindow(bool* p_open) {
    ImGui::SetNextWindowSizeConstraints(ImVec2(635, 0), ImVec2(635, FLT_MAX));
    if (!ImGui::Begin("Timer Window", p_open))
    {
        ImGui::End();
        return;
    }
    int32_t i = 0;
    for (auto tptr : TimerManager::getTimers()) {
        auto& t = *tptr;
        ImGui::PushStyleColor(ImGuiCol_Header, hexToCol(t.getColor()));
        auto id = t.getDecriptor() + (t.getSamples().size() > 0 ? 
            ("\t" + std::to_string(std::floor(*t.getSamples().rbegin() * 100.f) * 0.01f) + "ms") : std::string(""));
        if (ImGui::CollapsingHeader((id + "###" + (t.getDecriptor() + std::to_string(i++))).c_str())) {
            ImGui::PushID(i);
            static char buf1[512] = "";
            strcpy(buf1, t.getDecriptor().c_str());
            ImGui::PushItemWidth(-100);
            ImGui::InputText("Identifier", buf1, 64);
            if (t.getSourceLocation() != "") {
                strcpy(buf1, (t.getSourceLocation() + ":" + std::to_string(t.getLine())).c_str());
                ImGui::InputText("Location", buf1, 64);
            }
            ImGui::PopItemWidth();
            if (auto statsopt = t.getStats()) {
                auto stats = statsopt.value();
                static char buf[512] = "";
                sprintf(buf, "%.7g", std::floor(stats.min * 100.f) * 0.01f);
                ImGui::Text("min"); ImGui::SameLine();
                ImGui::Button(buf, ImVec2(66, 0)); ImGui::SameLine();
                ImGui::Text("max"); ImGui::SameLine();
                sprintf(buf, "%.7g", std::floor(stats.max * 100.f) * 0.01f);
                ImGui::Button(buf, ImVec2(66, 0)); ImGui::SameLine();
                ImGui::Text("avg"); ImGui::SameLine();
                sprintf(buf, "%.7g", std::floor(stats.avg * 100.f) * 0.01f);
                ImGui::Button(buf, ImVec2(66, 0)); ImGui::SameLine();
                ImGui::Text("med"); ImGui::SameLine();
                sprintf(buf, "%.7g", std::floor(stats.median * 100.f) * 0.01f);
                ImGui::Button(buf, ImVec2(66, 0)); ImGui::SameLine();
                ImGui::Text("dev"); ImGui::SameLine();
                sprintf(buf, "%.7g", std::floor(stats.stddev * 100.f) * 0.01f);
                ImGui::Button(buf, ImVec2(66, 0));

                ImGui::PushItemWidth(-100);
                ImGui::PlotLines("Timings", t.getSamples().data(), t.getSamples().size(), 0.f, NULL, stats.min, stats.max, ImVec2(0, 80));
                auto bmax = 0.f;
                for (auto& b : stats.histogram)
                    bmax = std::max(b, bmax);
                ImGui::PlotHistogram("Histogram", stats.histogram.data(), stats.histogram.size(), 0, NULL, 0, bmax, ImVec2(0, 80));
                ImGui::PopItemWidth();
            }

            ImGui::PopID();
        }
        ImGui::PopStyleColor();
    }
    ImGui::End();
}
