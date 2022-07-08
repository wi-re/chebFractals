#include <glad/glad.h>
#include "glui.h"
#include <imgui/implot.h>
#include <fsVis/utility.h>
#include <fsVis/optimizers.h>

void GUI::RayWindow(bool *p_open) {
  // return;
  ImGui::SetNextWindowSizeConstraints(ImVec2(1024, 512), ImVec2(FLT_MAX, FLT_MAX));
  if (!ImGui::Begin("Cheb Window", p_open)) {
    ImGui::End();
    return;
  }
  static std::array<float, 1024> xValues;
  static std::array<float, 1024> yValues;
  std::vector<float> coefficients;
  static bool once = true;
  float xmin = -1.;
  float xmax = 1.;
  float fmin = 0.;
  float fmax = 1.;

  {
    for (int32_t i = 0; i < 1024; ++i) {
      scalar x = -1. + ((scalar)i) / 1023. * 2.;
      xValues[i] = (float)x;
      yValues[i] = (float)globalFunction(x);
    }
    fmin = *std::min_element(yValues.begin(), yValues.begin() + 1024);
    fmax = *std::max_element(yValues.begin(), yValues.begin() + 1024);

    for (const auto &coeff : globalFunction.funs[0].coeffs())
      coefficients.push_back(std::abs((float)coeff));

    once = false;
  }

  if (ImPlot::BeginPlot("Absolute plot", ImVec2(1024, 256))) {
    ImPlot::SetupAxes("x", "f(x)");
    ImPlot::PushStyleVar(ImPlotStyleVar_LineWeight, 2.5f);
    ImPlot::PlotLine("Function", xValues.data(), yValues.data(), 1024);
    ImPlot::PopStyleVar();
    ImPlot::PushStyleVar(ImPlotStyleVar_LineWeight, 2.5f);
    ImPlot::PopStyleVar();
    ImPlot::EndPlot();
  }
  if (ImPlot::BeginPlot("Coefficient plot", ImVec2(1024, 256))) {
    ImPlot::SetupAxes("x", "f(x)");
    ImPlot::SetupAxis(ImAxis_Y1, NULL, ImPlotAxisFlags_LogScale);
    ImPlot::PushStyleVar(ImPlotStyleVar_LineWeight, 2.5f);
    ImPlot::PlotLine("f", coefficients.data(), coefficients.size());
    ImPlot::PopStyleVar();
    ImPlot::PushStyleVar(ImPlotStyleVar_LineWeight, 2.5f);
    ImPlot::PopStyleVar();
    ImPlot::EndPlot();
  }
  ImGui::End();
}
