#include <glad/glad.h> 
#include "glui.h"
#include <tools/timer.h>
#include <imgui/implot.h>
#include <simulation/SPHRender.h>

void GUI::RayWindow(bool* p_open) {
    return;
    ImGui::SetNextWindowSizeConstraints(ImVec2(1024, 512), ImVec2(FLT_MAX, FLT_MAX));
    if (!ImGui::Begin("Cheb Window", p_open))
    {
        ImGui::End();
        return;
    }
    static std::array<float, 1024> xValues;
    static std::array<float, 1024> yValues;
    static std::vector<float> coefficients;
    static std::vector<float> coefficientsdx;
    static std::vector<float> coefficientsd2x;
    static bool once = true;
    static float xmin = -1.;
    static float xmax = 1.;
    static float fmin = 0.;
    static float fmax = 1.;

    if(once){
    for(int32_t i = 0; i < 1024; ++i){
        scalar x = -1. + ((scalar) i) / 1023. * 2.;
        xValues[i] = (float) x;
        yValues[i] = (float) globalFunction(x);
    }
    fmin = *std::min_element(yValues.begin(), yValues.begin() + 1024);
    fmax = *std::max_element(yValues.begin(), yValues.begin() + 1024);

    for(const auto& coeff : globalFunction.funs[0].coeffs())
    coefficients.push_back(std::abs((float)coeff));
    for(const auto& coeff : globalFunctionFirstDerivative.funs[0].coeffs())
    coefficientsdx.push_back(std::abs((float)coeff));
    for(const auto& coeff : globalFunctionSecondDerivative.funs[0].coeffs())
    coefficientsd2x.push_back(std::abs((float)coeff));

    once = false;
    }


        // ImPlot::SetNextAxisLimits(ImAxis_X1, xmin, xmax, ImGuiCond_Always);
        // ImPlot::SetNextAxisLimits(ImAxis_Y1, fmin, fmax, ImGuiCond_Always);

    if (ImPlot::BeginPlot("Absolute plot", ImVec2(1024,256))) {
        ImPlot::SetupAxes("x","f(x)");
        ImPlot::PushStyleVar(ImPlotStyleVar_LineWeight, 2.5f);
        ImPlot::PlotLine("Function", xValues.data(), yValues.data(), 1024);
        ImPlot::PopStyleVar();
        ImPlot::PushStyleVar(ImPlotStyleVar_LineWeight, 2.5f);
        // ImPlot::PlotLine("Explicit Base", angles.data(), distancesExplicit.data(), angles.size());
        // ImPlot::PlotLine("Explicit Fine", angles.data(), distancesExplicitFine.data(), angles.size());
        ImPlot::PopStyleVar();
        ImPlot::EndPlot();
    }
    if (ImPlot::BeginPlot("Coefficient plot",ImVec2(1024,256))) {
        ImPlot::SetupAxes("x","f(x)");
        ImPlot::SetupAxis(ImAxis_Y1, NULL, ImPlotAxisFlags_LogScale);
        ImPlot::PushStyleVar(ImPlotStyleVar_LineWeight, 2.5f);
        ImPlot::PlotLine("f", coefficients.data(), coefficients.size());
        ImPlot::PlotLine("f'", coefficientsdx.data(), coefficientsdx.size());
        ImPlot::PlotLine("f''", coefficientsd2x.data(), coefficientsd2x.size());
        ImPlot::PopStyleVar();
        ImPlot::PushStyleVar(ImPlotStyleVar_LineWeight, 2.5f);
        // ImPlot::PlotLine("Explicit Base", angles.data(), distancesExplicit.data(), angles.size());
        // ImPlot::PlotLine("Explicit Fine", angles.data(), distancesExplicitFine.data(), angles.size());
        ImPlot::PopStyleVar();
        ImPlot::EndPlot();
    }
    // std::vector<float> diffFine, diffBase;
    // diffFine.resize(distancesExplicit.size());
    // diffBase.resize(distancesExplicit.size());
    // float minDist = FLT_MAX;
    // float maxDist = -FLT_MAX;
    // for (auto i = 0; i < distancesExplicit.size(); ++i) {
    //     diffFine[i] = distancesImplicit[i] - distancesExplicitFine[i];
    //     diffBase[i] = distancesImplicit[i] - distancesExplicit[i];
    //     minDist = std::min(minDist, diffFine[i]);
    //     maxDist = std::max(maxDist, diffFine[i]);
    //     minDist = std::min(minDist, diffBase[i]);
    //     maxDist = std::max(maxDist, diffBase[i]);
    // }
    // ImPlot::SetNextPlotLimits(-degToRad(fov / 2.0), degToRad(fov / 2.0), minDist - 1e-5f, maxDist + 1e-5f, ImGuiCond_Always);
    // if (ImPlot::BeginPlot("Error Plot", "x", "f(x)", ImVec2(1900, 500))) {
    //     ImPlot::PushStyleVar(ImPlotStyleVar_LineWeight, 2.5f);
    //     ImPlot::PlotLine("Implicit - Fine", angles.data(), diffFine.data(), angles.size());
    //     ImPlot::PlotLine("Implicit - Base", angles.data(), diffBase.data(), angles.size());
    //     ImPlot::PopStyleVar();
    //     ImPlot::EndPlot();
    // }

    //for (auto tptr : TimerManager::getTimers()) {
    //    auto& t = *tptr;
    //    ImGui::PushStyleColor(ImGuiCol_Header, hexToCol(t.getColor()));
    //    auto id = t.getDecriptor() + (t.getSamples().size() > 0 ? 
    //        ("\t" + std::to_string(std::floorf(*t.getSamples().rbegin() * 100.f) * 0.01f) + "ms") : std::string(""));
    //    if (ImGui::CollapsingHeader((id + "###" + (t.getDecriptor() + std::to_string(i++))).c_str())) {
    //        ImGui::PushID(i);
    //        static char buf1[512] = "";
    //        strcpy_s(buf1, 512, t.getDecriptor().c_str());
    //        ImGui::PushItemWidth(-100);
    //        ImGui::InputText("Identifier", buf1, 64);
    //        if (t.getSourceLocation() != "") {
    //            strcpy_s(buf1, 512, (t.getSourceLocation() + ":" + std::to_string(t.getLine())).c_str());
    //            ImGui::InputText("Location", buf1, 64);
    //        }
    //        ImGui::PopItemWidth();
    //        if (auto statsopt = t.getStats()) {
    //            auto stats = statsopt.value();
    //            static char buf[512] = "";
    //            sprintf_s(buf, 512, "%.7g", std::floorf(stats.min * 100.f) * 0.01f);
    //            ImGui::Text("min"); ImGui::SameLine();
    //            ImGui::Button(buf, ImVec2(66, 0)); ImGui::SameLine();
    //            ImGui::Text("max"); ImGui::SameLine();
    //            sprintf_s(buf, 512, "%.7g", std::floorf(stats.max * 100.f) * 0.01f);
    //            ImGui::Button(buf, ImVec2(66, 0)); ImGui::SameLine();
    //            ImGui::Text("avg"); ImGui::SameLine();
    //            sprintf_s(buf, 512, "%.7g", std::floorf(stats.avg * 100.f) * 0.01f);
    //            ImGui::Button(buf, ImVec2(66, 0)); ImGui::SameLine();
    //            ImGui::Text("med"); ImGui::SameLine();
    //            sprintf_s(buf, 512, "%.7g", std::floorf(stats.median * 100.f) * 0.01f);
    //            ImGui::Button(buf, ImVec2(66, 0)); ImGui::SameLine();
    //            ImGui::Text("dev"); ImGui::SameLine();
    //            sprintf_s(buf, 512, "%.7g", std::floorf(stats.stddev * 100.f) * 0.01f);
    //            ImGui::Button(buf, ImVec2(66, 0));

    //            ImGui::PushItemWidth(-100);
    //            ImGui::PlotLines("Timings", t.getSamples().data(), t.getSamples().size(), 0.f, NULL, stats.min, stats.max, ImVec2(0, 80));
    //            auto bmax = 0.f;
    //            for (auto& b : stats.histogram)
    //                bmax = std::max(b, bmax);
    //            ImGui::PlotHistogram("Histogram", stats.histogram.data(), stats.histogram.size(), 0, NULL, 0, bmax, ImVec2(0, 80));
    //            ImGui::PopItemWidth();
    //        }

    //        ImGui::PopID();
    //    }
    //    ImGui::PopStyleColor();
    //}
    ImGui::End();
}
