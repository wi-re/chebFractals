// #define _CRT_SECURE_NO_WARNINGS
#define _SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING
#include <tools/ParameterManager.h>
#include "../imgui/imgui.h"
#include <glad/glad.h> 
#include <GLFW/glfw3.h>
#pragma warning(push)
#pragma warning(disable:4251; disable:4275)
#include <yaml-cpp/yaml.h>
#pragma warning(pop)
//#include <utility/identifier/resource_helper.h>
#include <glm/glm.hpp>

// struct float2{float x,y;};
// struct float3{float x,y,z;};
// struct float4{float x,y,z,w;};
// struct uint2{uint x,y;};
// struct uint3{uint x,y,z;};
// struct uint4{uint x,y,z,w;};
// struct int2{int x,y;};
// struct int3{int x,y,z;};
// struct int4{int x,y,z,w;};
// struct double2{double x,y;};
// struct double3{double x,y,z;};
// struct double4{double x,y,z,w;};


namespace YAML {
#define CONVERT_V1(ty, cty)
#define CONVERT_V2(ty, cty)\
	template<>\
	struct convert<ty ## 2> {\
		static Node encode(const ty ## 2 & rhs) {\
			Node node;\
			node.push_back(rhs.x());\
			node.push_back(rhs.y());\
			node.SetStyle(YAML::EmitterStyle::value::Flow);\
			return node;\
		}\
		static bool decode(const Node& node, ty ## 2& rhs) {\
			if (!node.IsSequence() || node.size() != 2) {\
				return false;\
			}\
			rhs.x() = node[0].as<cty>();\
			rhs.y() = node[1].as<cty>();\
			return true;\
		}\
	};
#define CONVERT_V3(ty, cty)\
	template<>\
	struct convert<ty ## 3> {\
		static Node encode(const ty ## 3 & rhs) {\
			Node node;\
			node.push_back(rhs.x());\
			node.push_back(rhs.y());\
			node.push_back(rhs.z());\
			node.SetStyle(YAML::EmitterStyle::value::Flow);\
			return node;\
		}\
		static bool decode(const Node& node, ty ## 3& rhs) {\
			if (!node.IsSequence() || node.size() != 3) {\
				return false;\
			}\
			rhs.x() = node[0].as<cty>();\
			rhs.y() = node[1].as<cty>();\
			rhs.z() = node[2].as<cty>();\
			return true;\
		}\
	};
#define CONVERT_V4(ty, cty)\
	template<>\
	struct convert<ty ## 4> {\
		static Node encode(const ty ## 4 & rhs) {\
			Node node;\
			node.push_back(rhs.x());\
			node.push_back(rhs.y());\
			node.push_back(rhs.z());\
			node.push_back(rhs.w());\
			node.SetStyle(YAML::EmitterStyle::value::Flow);\
			return node;\
		}\
		static bool decode(const Node& node, ty ## 4& rhs) {\
			if (!node.IsSequence() || node.size() != 4) {\
				return false;\
			}\
			rhs.x() = node[0].as<cty>();\
			rhs.y() = node[1].as<cty>();\
			rhs.z() = node[2].as<cty>();\
			rhs.w() = node[3].as<cty>();\
			return true;\
		}\
	};
	template<>
	struct convert<glm::mat4> {
		static Node encode(const glm::mat4& rhs) {
			Node node;
			node.push_back(rhs[0][0]);
			node.push_back(rhs[0][1]);
			node.push_back(rhs[0][2]);
			node.push_back(rhs[0][3]);
			node.push_back(rhs[1][0]);
			node.push_back(rhs[1][1]);
			node.push_back(rhs[1][2]);
			node.push_back(rhs[1][3]);
			node.push_back(rhs[2][0]);
			node.push_back(rhs[2][1]);
			node.push_back(rhs[2][2]);
			node.push_back(rhs[2][3]);
			node.push_back(rhs[3][0]);
			node.push_back(rhs[3][1]);
			node.push_back(rhs[3][2]);
			node.push_back(rhs[3][3]);
			node.SetStyle(YAML::EmitterStyle::value::Flow);
			return node;
		}
		static bool decode(const Node& node, glm::mat4& rhs) {
			if (!node.IsSequence() || node.size() != 16) {
				return false;
			}
			rhs[0][0] = node[0].as<float>();
			rhs[0][1] = node[1].as<float>();
			rhs[0][2] = node[2].as<float>();
			rhs[0][3] = node[3].as<float>();
			rhs[1][0] = node[0].as<float>();
			rhs[1][1] = node[1].as<float>();
			rhs[1][2] = node[2].as<float>();
			rhs[1][3] = node[3].as<float>();
			rhs[2][0] = node[0].as<float>();
			rhs[2][1] = node[1].as<float>();
			rhs[2][2] = node[2].as<float>();
			rhs[2][3] = node[3].as<float>();
			rhs[3][0] = node[0].as<float>();
			rhs[3][1] = node[1].as<float>();
			rhs[3][2] = node[2].as<float>();
			rhs[3][3] = node[3].as<float>();
			return true;
		}
	};
#define CONVERTER(ty, cty) CONVERT_V1(ty, cty) CONVERT_V2(ty, cty) CONVERT_V3(ty, cty) CONVERT_V4(ty, cty)
	CONVERTER(float, float);
	CONVERTER(double, double);
	CONVERTER(int, int32_t);
}

namespace detail {
	std::string to_string(float4 rhs) {
		return "[ " + std::to_string(rhs.x()) + ", " + std::to_string(rhs.y()) + ", " + std::to_string(rhs.z()) + ", " + std::to_string(rhs.w()) + "]";
	}
	std::string to_string(double4 rhs) {
		return "[ " + std::to_string(rhs.x()) + ", " + std::to_string(rhs.y()) + ", " + std::to_string(rhs.z()) + ", " + std::to_string(rhs.w()) + "]";
	}
	std::string to_string(int4 rhs) {
		return "[ " + std::to_string(rhs.x()) + ", " + std::to_string(rhs.y()) + ", " + std::to_string(rhs.z()) + ", " + std::to_string(rhs.w()) + "]";
	}

	std::string to_string(float3 rhs) {
		return "[ " + std::to_string(rhs.x()) + ", " + std::to_string(rhs.y()) + ", " + std::to_string(rhs.z()) + "]";
	}
	std::string to_string(double3 rhs) {
		return "[ " + std::to_string(rhs.x()) + ", " + std::to_string(rhs.y()) + ", " + std::to_string(rhs.z()) + "]";
	}
	std::string to_string(int3 rhs) {
		return "[ " + std::to_string(rhs.x()) + ", " + std::to_string(rhs.y()) + ", " + std::to_string(rhs.z()) + "]";
	}

	std::string to_string(float2 rhs) {
		return "[ " + std::to_string(rhs.x()) + ", " + std::to_string(rhs.y()) + "]";
	}
	std::string to_string(double2 rhs) {
		return "[ " + std::to_string(rhs.x()) + ", " + std::to_string(rhs.y()) + "]";
	}
	std::string to_string(int2 rhs) {
		return "[ " + std::to_string(rhs.x()) + ", " + std::to_string(rhs.y()) + "]";
	}
}


std::pair<std::string, std::string> split(std::string s) {
	if (s.find(".") == std::string::npos)
		return std::make_pair(std::string(""), s);
	return std::make_pair(s.substr(0, s.find(".")), s.substr(s.find(".") + 1));
}

bool ParameterManager::parameterExists(std::string s) {
	if (parameterList.find(s) == parameterList.end()) return false;
	return true;
}
bool ParameterManager::isAmbiguous(std::string s) {
	auto idc = qualifiedIdentifierList.count(s);
	if (idc == 0)
		return false;
	//throw std::invalid_argument("Parameter " + id + " does not exist");
	if (idc > 1)
		return true;
	//throw std::invalid_argument("Parameter " + id + " is ambiguous");
	return false;
	//auto range = qualifiedIdentifierList.equal_range(id);
	//auto qIdentifier = range.first->second;
	//if (idc > 1)
	//	for (auto i = range.first; i != range.second; ++i)
	//		if (ns == i->second.substr(0, i->second.find(".")))
	//			qIdentifier = i->second;
	//return qIdentifier;
}
std::string ParameterManager::resolveParameter(std::string s) {
	if (s.find(".") != std::string::npos) return s;
	auto [ns, id] = split(s);
	auto idc = qualifiedIdentifierList.count(id);
	if (idc == 0)
		throw std::invalid_argument("Parameter " + id + " does not exist");
	if (idc > 1 && ns == "")
		throw std::invalid_argument("Parameter " + id + " is ambiguous");
	auto range = qualifiedIdentifierList.equal_range(id);
	auto qIdentifier = range.first->second;
	if (idc > 1)
		for (auto i = range.first; i != range.second; ++i)
			if (ns == i->second.substr(0, i->second.find(".")))
				qIdentifier = i->second;
	return qIdentifier;
}
void ParameterManager::parseTree(YAML::Node root) {
	for (auto& p : parameterList) {
		auto id = p.second->identifier;
		auto ns = p.second->identifierNamespace;
		try {
			YAML::Node node;
			bool found = false;
			if (ns != "") {
				for (auto nnode : root) {
					//std::cout << nnode.first << " : " << nnode.first.as<std::string>() << std::endl;
					//std::cout << nnode.first.IsSequence() << std::endl;
					//std::cout << nnode.first.IsMap() << std::endl;
					if (nnode.first.as<std::string>() == ns)
						if (nnode.second[id]) {
							node = nnode.second[id];
							found = true;
						}
				}
				//if (root[ns]) {
				//	auto nns = root[ns];
				//	if (nns[id])
				//		node = nns[id];
				//}
			}
			else {
				node = root[id];
				found = true;
			}
			if (found && node) {
				if (!p.second->param.isVec) {
					if (p.second->param.valVec)
						p.second->param.val.value() = decoders[p.second->type](node);
					else
						p.second->param.val = decoders[p.second->type](node);
				}
				else {
					if (!node.IsSequence() && !node.IsNull())
						throw std::invalid_argument("Expected sequence for " + p.first);
					if (p.second->param.valVec)
						p.second->param.valVec.value() = decoders[p.second->type](node);
					else
						p.second->param.valVec = decoders[p.second->type](node);
				}
			}
		}
		catch (...) {
			std::cout << ns << " : " << id << std::endl;
			throw;
		}
	}
}

#define customVectorInternal(ty)\
	addUifunction(typeid(std::vector<ty>), [](Parameter& parameter) {\
		std::vector<ty>& vec = boost::any_cast<std::vector<ty>&>(parameter.param.valVec.value());\
		if (parameter.properties.hidden) return;\
		ImGui::PushID(parameter.identifier.c_str());\
		ImGui::Text("%s", parameter.identifier.c_str());\
		if (!parameter.properties.constant) {\
			ImGui::SameLine();\
			if (ImGui::Button("+"))\
				vec.push_back(ty());\
			ImGui::SameLine();\
			if (ImGui::Button("-"))\
				vec.pop_back();\
		}\
		int32_t i = 0;\
		for (auto& elem : vec) {\
			ImGui::Indent();\
			ImGui::PushID(i);\
			auto eParam = Parameter{ parameter.identifier + "[" + std::to_string(i++) + "]", parameter.identifierNamespace, elem, typeid(ty), parameter.properties };\
			ParameterManager::instance().uiFunctions[typeid(ty)](eParam);\
			ImGui::Unindent();\
			ImGui::PopID();\
			if (!parameter.properties.constant)\
				elem = (ty) eParam.param.val.value();\
		}\
		ImGui::PopID();\
		});


ParameterManager::ParameterManager() {
#define SIMPLE_PARSER(ty) \
		decoders[typeid(ty)] = [](const YAML::Node& node) {return detail::iAny(node.as<ty>()); };\
		encoders[typeid(ty)] = [](const detail::iAny& any) {return YAML::convert<ty>::encode(boost::any_cast<ty>(any)); };
#define VECTOR_PARSER(ty)\
		decoders[typeid(std::vector<ty>)] = [](const YAML::Node& node) {\
			std::vector<ty> vec;\
			for (auto e : node)\
				vec.push_back(e.as<ty>());\
			return detail::iAny(vec); \
		};\
		encoders[typeid(std::vector<ty>)] = [](const detail::iAny& any){return YAML::convert<std::vector<ty>>::encode(boost::any_cast<std::vector<ty>>(any));};
	SIMPLE_PARSER(int32_t);
	SIMPLE_PARSER(uint32_t);
	SIMPLE_PARSER(float);
	SIMPLE_PARSER(double);
	SIMPLE_PARSER(bool);
	SIMPLE_PARSER(std::string);
	SIMPLE_PARSER(float2);
	SIMPLE_PARSER(float3);
	SIMPLE_PARSER(float4);
	SIMPLE_PARSER(double2);
	SIMPLE_PARSER(double3);
	SIMPLE_PARSER(double4);
	SIMPLE_PARSER(int2);
	SIMPLE_PARSER(int3);
	SIMPLE_PARSER(int4);
	SIMPLE_PARSER(glm::mat4);
	VECTOR_PARSER(int32_t);
	VECTOR_PARSER(uint32_t);
	VECTOR_PARSER(float);
	VECTOR_PARSER(double);
	VECTOR_PARSER(std::string);
	VECTOR_PARSER(float2);
	VECTOR_PARSER(float3);
	VECTOR_PARSER(float4);
	VECTOR_PARSER(double2);
	VECTOR_PARSER(double3);
	VECTOR_PARSER(double4);
	VECTOR_PARSER(int2);
	VECTOR_PARSER(int3);
	VECTOR_PARSER(int4);

	decoders[typeid(char const*)] = [](const YAML::Node& node) {
		auto str = node.as<std::string>();
		char* pt = new char[str.size() + 1];
		memcpy(pt, str.c_str(), str.size() + 1);
		char const* cpt = pt;
		return detail::iAny(cpt); };
	encoders[typeid(char const*)] = [](const detail::iAny& node) {
		char const* cpt = boost::any_cast<char const*>(node);
		return YAML::convert < std::string>::encode(std::string(cpt));
	};
	addUifunction(typeid(bool), [](Parameter& param) {
		bool& var = boost::any_cast<bool&>(param.param.val.value());
		//std::cout << param.identifier << std::endl;
		if (param.properties.hidden)return;
		if (param.properties.constant) {
			int32_t ib = var ? 1 : 0;
			auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
			ImGui::SliderInt(param.identifier.c_str(), &ib, 0, 1);
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
			if (param.properties.description != "" && ImGui::IsItemHovered())
				ImGui::SetTooltip("%s", param.properties.description.c_str());
			return;
		}
		int32_t ib = var ? 1 : 0;
		int32_t ibb = ib;
		ImGui::SliderInt(param.identifier.c_str(), &ib, 0, 1);
		if (param.properties.description != "" && ImGui::IsItemHovered())
			ImGui::SetTooltip("%s", param.properties.description.c_str());
		if (ib != ibb)
			var = ib == 0 ? false : true;
		return;
		});
	addUifunction(typeid(std::string), [](Parameter& param) {
		std::string& var = boost::any_cast<std::string&>(param.param.val.value());
		if (param.properties.hidden)return;
		if (param.properties.constant) {
			auto vcp = var;
			auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
			static char buf1[256] = "";
			#ifdef WIN32
			strcpy_s(buf1, 256, var.c_str());
			#else
			strcpy(buf1,var.c_str());
			#endif
			ImGui::InputText(param.identifier.c_str(), buf1, 64);
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
			if (param.properties.description != "" && ImGui::IsItemHovered())
				ImGui::SetTooltip("%s", param.properties.description.c_str());
			return;
		}
		static char buf1[256] = "";
			#ifdef WIN32
			strcpy_s(buf1, 256, var.c_str());
			#else
			strcpy(buf1,var.c_str());
			#endif
		ImGui::InputText(param.identifier.c_str(), buf1, 64);
		var = buf1;

		if (param.properties.description != "" && ImGui::IsItemHovered())
			ImGui::SetTooltip("%s", param.properties.description.c_str());
		if (param.properties.presets.size() != 0) {
			std::vector<std::string> presets;
			for (auto& pr : param.properties.presets)
				presets.push_back(boost::any_cast<std::string>(pr));
			ImGui::SameLine();
			if (ImGui::BeginCombo((param.identifier + " presets").c_str(), var.c_str())) {
				for (int n = 0; n < presets.size(); n++) {
					bool is_selected = (var == presets[n]);
					if (ImGui::Selectable(presets[n].c_str(), is_selected))
						var = presets[n];
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
		}
		});
	addUifunction(typeid(char const*), [](Parameter& param) {
		char const* var = boost::any_cast<char const*>(param.param.val.value());
		if (param.properties.hidden)return;
		auto vcp = var;
		auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
		ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
		static char buf1[256] = "";
			#ifdef WIN32
			strcpy_s(buf1, 256, var);
			#else
			strcpy(buf1,var);
			#endif
		ImGui::InputText(param.identifier.c_str(), buf1, 64);
		ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
		if (param.properties.description != "" && ImGui::IsItemHovered())
			ImGui::SetTooltip("%s", param.properties.description.c_str());
		return;
		});
	addUifunction(typeid(int32_t), [](Parameter& param) {
		int32_t& var = boost::any_cast<int32_t&>(param.param.val.value());
		if (param.properties.hidden)return;
		if (param.properties.constant) {
			auto vcp = var;
			auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
			ImGui::DragInt(param.identifier.c_str(), &vcp, 0, vcp, vcp);
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
			if (param.properties.description != "" && ImGui::IsItemHovered())
				ImGui::SetTooltip("%s", param.properties.description.c_str());
			return;
		}
		if (param.properties.range)
			ImGui::SliderInt(param.identifier.c_str(), &var, param.properties.range.value().min, param.properties.range.value().max);
		else ImGui::DragInt(param.identifier.c_str(), &var);
		if (param.properties.description != "" && ImGui::IsItemHovered())
			ImGui::SetTooltip("%s", param.properties.description.c_str());
		if (param.properties.presets.size() != 0) {
			std::vector<int32_t> presets;
			for (auto& pr : param.properties.presets)
				presets.push_back(boost::any_cast<int32_t>(pr));
			ImGui::SameLine();
			if (ImGui::BeginCombo((param.identifier + " presets").c_str(), std::to_string(var).c_str())) {
				for (int n = 0; n < presets.size(); n++) {
					bool is_selected = (var == presets[n]);
					if (ImGui::Selectable(std::to_string(presets[n]).c_str(), is_selected))
						var = presets[n];
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
		}
		});
	addUifunction(typeid(uint32_t), [](Parameter& param) {
		uint32_t& var = boost::any_cast<uint32_t&>(param.param.val.value());
		int32_t vari = (int32_t)var;
		int32_t varii = (int32_t)var;
		if (param.properties.hidden)return;
		if (param.properties.constant) {
			auto vcp = vari;
			auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
			ImGui::DragInt(param.identifier.c_str(), &vcp, 0, vcp, vcp);
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
			//ImGui::Text((param.identifier + ": " + std::to_string(var)).c_str());
			if (param.properties.description != "" && ImGui::IsItemHovered())
				ImGui::SetTooltip("%s", param.properties.description.c_str());
			return;
		}
		if (param.properties.range)
			ImGui::SliderInt(param.identifier.c_str(), &vari, (int32_t)(uint32_t)param.properties.range.value().min, (int32_t)(uint32_t)param.properties.range.value().max);
		else ImGui::DragInt(param.identifier.c_str(), &vari, 1, 0, INT_MAX);
		if (param.properties.description != "" && ImGui::IsItemHovered())
			ImGui::SetTooltip("%s", param.properties.description.c_str());
		if (varii != vari) var = (uint32_t)vari;
		if (param.properties.presets.size() != 0) {
			std::vector<uint32_t> presets;
			for (auto& pr : param.properties.presets)
				presets.push_back(boost::any_cast<uint32_t>(pr));
			ImGui::SameLine();
			if (ImGui::BeginCombo((param.identifier + " presets").c_str(), std::to_string(var).c_str())) {
				for (int n = 0; n < presets.size(); n++) {
					bool is_selected = (var == presets[n]);
					if (ImGui::Selectable(std::to_string(presets[n]).c_str(), is_selected))
						var = presets[n];
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
		}
		});
	addUifunction(typeid(std::size_t), [](Parameter& param) {
		std::size_t& var = boost::any_cast<std::size_t&>(param.param.val.value());
		int32_t vari = (int32_t)var;
		int32_t varii = (int32_t)var;
		if (param.properties.hidden)return;
		if (param.properties.constant) {
			auto vcp = vari;
			auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
			ImGui::DragInt(param.identifier.c_str(), &vcp, 0, vcp, vcp);
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
			//ImGui::Text((param.identifier + ": " + std::to_string(var)).c_str());
			if (param.properties.description != "" && ImGui::IsItemHovered())
				ImGui::SetTooltip("%s", param.properties.description.c_str());
			return;
		}
		if (param.properties.range)
			ImGui::SliderInt(param.identifier.c_str(), &vari, (int32_t)(std::size_t)param.properties.range.value().min, (int32_t)(std::size_t)param.properties.range.value().max);
		else ImGui::DragInt(param.identifier.c_str(), &vari, 1, 0, INT_MAX);
		if (param.properties.description != "" && ImGui::IsItemHovered())
			ImGui::SetTooltip("%s", param.properties.description.c_str());
		if (varii != vari) var = (std::size_t)vari;
		if (param.properties.presets.size() != 0) {
			std::vector<std::size_t> presets;
			for (auto& pr : param.properties.presets)
				presets.push_back(boost::any_cast<std::size_t>(pr));
			ImGui::SameLine();
			if (ImGui::BeginCombo((param.identifier + " presets").c_str(), std::to_string(var).c_str())) {
				for (int n = 0; n < presets.size(); n++) {
					bool is_selected = (var == presets[n]);
					if (ImGui::Selectable(std::to_string(presets[n]).c_str(), is_selected))
						var = presets[n];
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
		}
		});
	addUifunction(typeid(float), [](Parameter& param) {
		float& var = boost::any_cast<float&>(param.param.val.value());
		if (param.properties.hidden)return;
		if (param.properties.constant) {
			auto vcp = var;
			auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
			ImGui::DragFloat(param.identifier.c_str(), &vcp, 0, vcp, vcp, "%.6f");
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
			if (param.properties.description != "" && ImGui::IsItemHovered())
				ImGui::SetTooltip("%s", param.properties.description.c_str());
			return;
		}
		if (param.properties.range)
			ImGui::SliderFloat(param.identifier.c_str(), &var, param.properties.range.value().min, param.properties.range.value().max, "%.5g");
		else ImGui::DragFloat(param.identifier.c_str(), &var, var * 0.01f, -FLT_MAX, FLT_MAX, "%.6f");
		if (param.properties.description != "" && ImGui::IsItemHovered())
			ImGui::SetTooltip("%s", param.properties.description.c_str());
		if (param.properties.presets.size() != 0) {
			std::vector<float> presets;
			for (auto& pr : param.properties.presets)
				presets.push_back(boost::any_cast<float>(pr));
			ImGui::SameLine();
			if (ImGui::BeginCombo((param.identifier + " presets").c_str(), std::to_string(var).c_str())) {
				for (int n = 0; n < presets.size(); n++) {
					bool is_selected = (var == presets[n]);
					if (ImGui::Selectable(std::to_string(presets[n]).c_str(), is_selected))
						var = presets[n];
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
		}
		});
	addUifunction(typeid(double), [](Parameter& param) {
		double& var = boost::any_cast<double&>(param.param.val.value());
		float vard = (float)var;
		float vardd = vard;
		if (param.properties.hidden)return;
		if (param.properties.constant) {
			auto vcp = vard;
			auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
			ImGui::DragFloat(param.identifier.c_str(), &vcp, 0, vcp, vcp, "%.6f");
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
			if (param.properties.description != "" && ImGui::IsItemHovered())
				ImGui::SetTooltip("%s", param.properties.description.c_str());
			return;
		}
		if (param.properties.range)
			ImGui::SliderFloat(param.identifier.c_str(), &vard, (float)(double)param.properties.range.value().min, (float)(double)param.properties.range.value().max);
		else ImGui::DragFloat(param.identifier.c_str(), &vard, 0, 0.f, 0.f, "%.6f");
		if (param.properties.description != "" && ImGui::IsItemHovered())
			ImGui::SetTooltip("%s", param.properties.description.c_str());
		if (vard != vardd)
			var = vard;
		if (param.properties.presets.size() != 0) {
			std::vector<double> presets;
			for (auto& pr : param.properties.presets)
				presets.push_back(boost::any_cast<double>(pr));
			ImGui::SameLine();
			if (ImGui::BeginCombo((param.identifier + " presets").c_str(), std::to_string(var).c_str())) {
				for (int n = 0; n < presets.size(); n++) {
					bool is_selected = (var == presets[n]);
					if (ImGui::Selectable(std::to_string(presets[n]).c_str(), is_selected))
						var = presets[n];
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
		}
		});

	addUifunction(typeid(float4), [](Parameter& param) {
		float4& var = boost::any_cast<float4&>(param.param.val.value());
		if (param.properties.hidden)return;
		if (param.properties.constant) {
			auto vcp = var;
			auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
			ImGui::DragFloat4(param.identifier.c_str(), &vcp.x(), 0, 0, 0);
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
			if (param.properties.description != "" && ImGui::IsItemHovered())
				ImGui::SetTooltip("%s", param.properties.description.c_str());
			return;
		}
		if (param.properties.range)
			ImGui::SliderFloat4(param.identifier.c_str(), &var.x(), param.properties.range.value().min, param.properties.range.value().max);
		else ImGui::DragFloat4(param.identifier.c_str(), &var.x(), 0.01f);
		if (param.properties.description != "" && ImGui::IsItemHovered())
			ImGui::SetTooltip("%s", param.properties.description.c_str());
		if (param.properties.presets.size() != 0) {
			std::vector<float4> presets;
			for (auto& pr : param.properties.presets)
				presets.push_back(boost::any_cast<float4>(pr));
			ImGui::SameLine();
			if (ImGui::BeginCombo((param.identifier + " presets").c_str(), detail::to_string(var).c_str())) {
				for (int n = 0; n < presets.size(); n++) {
					bool is_selected = (detail::to_string(var) == detail::to_string(presets[n]));
					if (ImGui::Selectable(detail::to_string(presets[n]).c_str(), is_selected))
						var = presets[n];
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
		}
		});
	addUifunction(typeid(double4), [](Parameter& param) {
		double4& vard = boost::any_cast<double4&>(param.param.val.value());
		float4 var{ (float)vard.x(), (float)vard.y(), (float)vard.z(),(float)vard.w() };
		float4 varb = var;
		if (param.properties.hidden)return;
		if (param.properties.constant) {
			auto vcp = var;
			auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
			ImGui::DragFloat4(param.identifier.c_str(), &vcp.x(), 0, 0, 0);
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
			if (param.properties.description != "" && ImGui::IsItemHovered())
				ImGui::SetTooltip("%s", param.properties.description.c_str());
			return;
		}
		if (param.properties.range)
			ImGui::SliderFloat4(param.identifier.c_str(), &var.x(), (float)(double)param.properties.range.value().min, (float)(double)param.properties.range.value().max);
		else ImGui::DragFloat4(param.identifier.c_str(), &var.x(), 0.01f);
		if (param.properties.description != "" && ImGui::IsItemHovered())
			ImGui::SetTooltip("%s", param.properties.description.c_str());
		if (varb.x() != var.x() || varb.y() != var.y() || varb.z() != var.z() || varb.w() != var.w()) vard = double4{ (double)var.x(), (double)var.y(),(double)var.z(), (double)var.w() };
		if (param.properties.presets.size() != 0) {
			std::vector<double4> presets;
			for (auto& pr : param.properties.presets)
				presets.push_back(boost::any_cast<double4>(pr));
			ImGui::SameLine();
			if (ImGui::BeginCombo((param.identifier + " presets").c_str(), detail::to_string(var).c_str())) {
				for (int n = 0; n < presets.size(); n++) {
					bool is_selected = (detail::to_string(var) == detail::to_string(presets[n]));
					if (ImGui::Selectable(detail::to_string(presets[n]).c_str(), is_selected))
						vard = double4{ presets[n].x(), presets[n].y(), presets[n].z(), presets[n].w() };
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
		}
		});

	addUifunction(typeid(int4), [](Parameter& param) {
		int4& var = boost::any_cast<int4&>(param.param.val.value());
		if (param.properties.hidden)return;
		if (param.properties.constant) {
			auto vcp = var;
			auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
			ImGui::DragInt4(param.identifier.c_str(), &vcp.x(), 0, 0, 0);
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
			if (param.properties.description != "" && ImGui::IsItemHovered())
				ImGui::SetTooltip("%s", param.properties.description.c_str());
			return;
		}
		if (param.properties.range)
			ImGui::SliderInt4(param.identifier.c_str(), &var.x(), param.properties.range.value().min, param.properties.range.value().max);
		else ImGui::DragInt4(param.identifier.c_str(), &var.x(), 1.f);
		if (param.properties.description != "" && ImGui::IsItemHovered())
			ImGui::SetTooltip("%s", param.properties.description.c_str());
		if (param.properties.presets.size() != 0) {
			std::vector<int4> presets;
			for (auto& pr : param.properties.presets)
				presets.push_back(boost::any_cast<int4>(pr));
			ImGui::SameLine();
			if (ImGui::BeginCombo((param.identifier + " presets").c_str(), detail::to_string(var).c_str())) {
				for (int n = 0; n < presets.size(); n++) {
					bool is_selected = (detail::to_string(var) == detail::to_string(presets[n]));
					if (ImGui::Selectable(detail::to_string(presets[n]).c_str(), is_selected))
						var = presets[n];
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
		}
		});
	addUifunction(typeid(float3), [](Parameter& param) {
		float3& var = boost::any_cast<float3&>(param.param.val.value());
		if (param.properties.hidden)return;
		if (param.properties.constant) {
			auto vcp = var;
			auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
			ImGui::DragFloat3(param.identifier.c_str(), &vcp.x(), 0, 0, 0);
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
			if (param.properties.description != "" && ImGui::IsItemHovered())
				ImGui::SetTooltip("%s", param.properties.description.c_str());
			return;
		}
		if (param.properties.range)
			ImGui::SliderFloat3(param.identifier.c_str(), &var.x(), param.properties.range.value().min, param.properties.range.value().max);
		else ImGui::DragFloat3(param.identifier.c_str(), &var.x(), 0.01f);
		if (param.properties.description != "" && ImGui::IsItemHovered())
			ImGui::SetTooltip("%s", param.properties.description.c_str());
		if (param.properties.presets.size() != 0) {
			std::vector<float3> presets;
			for (auto& pr : param.properties.presets)
				presets.push_back(boost::any_cast<float3>(pr));
			ImGui::SameLine();
			if (ImGui::BeginCombo((param.identifier + " presets").c_str(), detail::to_string(var).c_str())) {
				for (int n = 0; n < presets.size(); n++) {
					bool is_selected = (detail::to_string(var) == detail::to_string(presets[n]));
					if (ImGui::Selectable(detail::to_string(presets[n]).c_str(), is_selected))
						var = presets[n];
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
		}
		});
	addUifunction(typeid(double3), [](Parameter& param) {
		double3& vard = boost::any_cast<double3&>(param.param.val.value());
		float3 var{ (float)vard.x(), (float)vard.y(), (float)vard.z() };
		float3 varb = var;
		if (param.properties.hidden)return;
		if (param.properties.constant) {
			auto vcp = var;
			auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
			ImGui::DragFloat3(param.identifier.c_str(), &vcp.x(), 0, 0, 0);
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
			if (param.properties.description != "" && ImGui::IsItemHovered())
				ImGui::SetTooltip("%s", param.properties.description.c_str());
			return;
		}
		if (param.properties.range)
			ImGui::SliderFloat3(param.identifier.c_str(), &var.x(), (float)(double)param.properties.range.value().min, (float)(double)param.properties.range.value().max);
		else ImGui::DragFloat3(param.identifier.c_str(), &var.x(), 0.01f);
		if (param.properties.description != "" && ImGui::IsItemHovered())
			ImGui::SetTooltip("%s", param.properties.description.c_str());
		if (varb.x() != var.x() || varb.y() != var.y() || varb.z() != var.z()) vard = double3{ (double)var.x(), (double)var.y(),(double)var.z() };
		if (param.properties.presets.size() != 0) {
			std::vector<double3> presets;
			for (auto& pr : param.properties.presets)
				presets.push_back(boost::any_cast<double3>(pr));
			ImGui::SameLine();
			if (ImGui::BeginCombo((param.identifier + " presets").c_str(), detail::to_string(var).c_str())) {
				for (int n = 0; n < presets.size(); n++) {
					bool is_selected = (detail::to_string(var) == detail::to_string(presets[n]));
					if (ImGui::Selectable(detail::to_string(presets[n]).c_str(), is_selected))
						vard = double3{ presets[n].x(), presets[n].y(), presets[n].z() };
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
		}
		});

	addUifunction(typeid(int3), [](Parameter& param) {
		int3& var = boost::any_cast<int3&>(param.param.val.value());
		if (param.properties.hidden)return;
		if (param.properties.constant) {
			auto vcp = var;
			auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
			ImGui::DragInt3(param.identifier.c_str(), &vcp.x(), 0, 0, 0);
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
			if (param.properties.description != "" && ImGui::IsItemHovered())
				ImGui::SetTooltip("%s", param.properties.description.c_str());
			return;
		}
		if (param.properties.range)
			ImGui::SliderInt3(param.identifier.c_str(), &var.x(), param.properties.range.value().min, param.properties.range.value().max);
		else ImGui::DragInt3(param.identifier.c_str(), &var.x(), 1.f);
		if (param.properties.description != "" && ImGui::IsItemHovered())
			ImGui::SetTooltip("%s", param.properties.description.c_str());
		if (param.properties.presets.size() != 0) {
			std::vector<int3> presets;
			for (auto& pr : param.properties.presets)
				presets.push_back(boost::any_cast<int3>(pr));
			ImGui::SameLine();
			if (ImGui::BeginCombo((param.identifier + " presets").c_str(), detail::to_string(var).c_str())) {
				for (int n = 0; n < presets.size(); n++) {
					bool is_selected = (detail::to_string(var) == detail::to_string(presets[n]));
					if (ImGui::Selectable(detail::to_string(presets[n]).c_str(), is_selected))
						var = presets[n];
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
		}
		});
	addUifunction(typeid(float2), [](Parameter& param) {
		float2& var = boost::any_cast<float2&>(param.param.val.value());
		if (param.properties.hidden)return;
		if (param.properties.constant) {
			auto vcp = var;
			auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
			ImGui::DragFloat2(param.identifier.c_str(), &vcp.x(), 0, 0, 0);
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
			if (param.properties.description != "" && ImGui::IsItemHovered())
				ImGui::SetTooltip("%s", param.properties.description.c_str());
			return;
		}
		if (param.properties.range)
			ImGui::SliderFloat2(param.identifier.c_str(), &var.x(), param.properties.range.value().min, param.properties.range.value().max);
		else ImGui::DragFloat2(param.identifier.c_str(), &var.x(), 0.01f);
		if (param.properties.description != "" && ImGui::IsItemHovered())
			ImGui::SetTooltip("%s", param.properties.description.c_str());
		if (param.properties.presets.size() != 0) {
			std::vector<float2> presets;
			for (auto& pr : param.properties.presets)
				presets.push_back(boost::any_cast<float2>(pr));
			ImGui::SameLine();
			if (ImGui::BeginCombo((param.identifier + " presets").c_str(), detail::to_string(var).c_str())) {
				for (int n = 0; n < presets.size(); n++) {
					bool is_selected = (detail::to_string(var) == detail::to_string(presets[n]));
					if (ImGui::Selectable(detail::to_string(presets[n]).c_str(), is_selected))
						var = presets[n];
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
		}
		});
	addUifunction(typeid(double2), [](Parameter& param) {
		double2& vard = boost::any_cast<double2&>(param.param.val.value());
		float2 var{ (float)vard.x(), (float)vard.y() };
		float2 varb = var;
		if (param.properties.hidden)return;
		if (param.properties.constant) {
			auto vcp = var;
			auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
			ImGui::DragFloat2(param.identifier.c_str(), &vcp.x(), 0, 0, 0);
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
			if (param.properties.description != "" && ImGui::IsItemHovered())
				ImGui::SetTooltip("%s", param.properties.description.c_str());
			return;
		}
		if (param.properties.range)
			ImGui::SliderFloat2(param.identifier.c_str(), &var.x(), (float)(double)param.properties.range.value().min, (float)(double)param.properties.range.value().max);
		else ImGui::DragFloat2(param.identifier.c_str(), &var.x(), 0.01f);
		if (param.properties.description != "" && ImGui::IsItemHovered())
			ImGui::SetTooltip("%s", param.properties.description.c_str());
		if (varb.x() != var.x() || varb.y() != var.y()) vard = double2{ (double)var.x(), (double)var.y() };
		if (param.properties.presets.size() != 0) {
			std::vector<double2> presets;
			for (auto& pr : param.properties.presets)
				presets.push_back(boost::any_cast<double2>(pr));
			ImGui::SameLine();
			if (ImGui::BeginCombo((param.identifier + " presets").c_str(), detail::to_string(var).c_str())) {
				for (int n = 0; n < presets.size(); n++) {
					bool is_selected = (detail::to_string(var) == detail::to_string(presets[n]));
					if (ImGui::Selectable(detail::to_string(presets[n]).c_str(), is_selected))
						vard = double2{ presets[n].x(), presets[n].y() };
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
		}
		});

	addUifunction(typeid(int2), [](Parameter& param) {
		int2& var = boost::any_cast<int2&>(param.param.val.value());
		if (param.properties.hidden)return;
		if (param.properties.constant) {
			auto vcp = var;
			auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
			ImGui::DragInt2(param.identifier.c_str(), &vcp.x(), 0, 0, 0);
			ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
			if (param.properties.description != "" && ImGui::IsItemHovered())
				ImGui::SetTooltip("%s", param.properties.description.c_str());
			return;
		}
		if (param.properties.range)
			ImGui::SliderInt2(param.identifier.c_str(), &var.x(), param.properties.range.value().min, param.properties.range.value().max);
		else ImGui::DragInt2(param.identifier.c_str(), &var.x(), 1.f);
		if (param.properties.description != "" && ImGui::IsItemHovered())
			ImGui::SetTooltip("%s", param.properties.description.c_str());
		if (param.properties.presets.size() != 0) {
			std::vector<int2> presets;
			for (auto& pr : param.properties.presets)
				presets.push_back(boost::any_cast<int2>(pr));
			ImGui::SameLine();
			if (ImGui::BeginCombo((param.identifier + " presets").c_str(), detail::to_string(var).c_str())) {
				for (int n = 0; n < presets.size(); n++) {
					bool is_selected = (detail::to_string(var) == detail::to_string(presets[n]));
					if (ImGui::Selectable(detail::to_string(presets[n]).c_str(), is_selected))
						var = presets[n];
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
		}
		});
	customVectorInternal(int32_t);
	customVectorInternal(float);
	/*	addUifunction(typeid(std::vector<int32_t>), [](Parameter& param) {
			std::vector<int32_t>& var = boost::any_cast<std::vector<int32_t>&>(param.param.valVec.value());
			if (param.properties.hidden)return;
			if (param.properties.constant)
			{
				ImGui::Text(param.identifier.c_str());
				int32_t i = 0;
				for (auto e : var) {
					ImGui::DragInt(("##" + param.identifier + std::to_string(i++)).c_str(), &e);
				}
			}
			ImGui::Text(param.identifier.c_str());
			ImGui::SameLine();
			if (ImGui::Button("+"))
				var.push_back(int32_t());
			ImGui::SameLine();
			if (ImGui::Button("-"))
				var.pop_back();
			int32_t i = 0;
			for (auto& e : var) {
				ImGui::DragInt(("##" + param.identifier + std::to_string(i++)).c_str(), &e);
			}

			});
		addUifunction(typeid(std::vector<float>), [](Parameter& param) {
			std::vector<float>& var = boost::any_cast<std::vector<float>&>(param.param.valVec.value());
			if (param.properties.hidden)return;
			if (param.properties.constant)
			{
				ImGui::Text(param.identifier.c_str());
				int32_t i = 0;
				for (auto e : var) {
					ImGui::DragFloat(("##" + param.identifier + std::to_string(i++)).c_str(), &e,e * 0.01f);
				}
			}
			ImGui::Text(param.identifier.c_str());
			ImGui::SameLine();
			if (ImGui::Button("+"))
				var.push_back(float());
			ImGui::SameLine();
			if (ImGui::Button("-"))
				var.pop_back();
			int32_t i = 0;
			for (auto& e : var) {
				ImGui::DragFloat(("##" + param.identifier + std::to_string(i++)).c_str(), &e, fabsf(e) > 1e-5f ? fabsf(e) * 0.01f : 1e-5f, -FLT_MAX, FLT_MAX);
			}

			});*/

}
void ParameterManager::addDecoder(std::type_index ty, std::function<detail::iAny(const YAML::Node&)> fn) { decoders[ty] = fn; }
void ParameterManager::addEncoder(std::type_index ty, std::function<YAML::Node(const detail::iAny&)> fn) { encoders[ty] = fn; }
void ParameterManager::addUifunction(std::type_index ty, std::function<void(Parameter&)> fn) { uiFunctions[ty] = fn; }
ParameterManager& ParameterManager::instance() {
	static ParameterManager inst;
	static bool once = true;
	return inst;
}
Parameter& ParameterManager::getParameter(std::string identifier) {
	auto qid = resolveParameter(identifier);
	if (!parameterExists(qid)) throw std::invalid_argument("Parameter " + identifier + " does not exist");
	return *parameterList[qid];
}
detail::varAny& ParameterManager::get(std::string identifier) {
	return getParameter(identifier).param;
}
void ParameterManager::load(std::string filename) {
	namespace fs = std::filesystem;
	if (!fs::exists(filename))
		throw std::invalid_argument("File does not exist");

	auto config = YAML::LoadFile(fs::path(filename).string());
	parseTree(config);
}
void ParameterManager::loadDirect(std::string yaml) {
	parseTree(YAML::Load(yaml));
}
void ParameterManager::loadTree(YAML::Node yaml) {
	parseTree(yaml);
}
YAML::Node ParameterManager::buildTree() {
	YAML::Node root;
	int32_t i = 0;
	for (auto p : parameterList) {
		auto id = p.second->identifier;
		auto ns = p.second->identifierNamespace;
		if (ns != "")
			root[ns][id] = encoders[p.second->type](p.second->param.isVec ? p.second->param.valVec.value() : p.second->param.val.value());
		else
			root[id] = encoders[p.second->type](p.second->param.isVec ? p.second->param.valVec.value() : p.second->param.val.value());
		//if (i++ > 50)break;
	}
	return root;
}
void ParameterManager::buildImguiWindow(bool* p_open) {
	{
		if (!ImGui::Begin("Parameter Manager", p_open))
		{
			ImGui::End();
			return;
		}
		std::map<std::string, std::vector<std::pair<std::string, Parameter*>>> parameters;
		for (auto p : parameterList) {
			auto [ns, id] = split(p.first);
			parameters[ns].push_back(std::make_pair(id, p.second));
		}
		for (auto param : parameters) {
			//std::cout << param.first << std::endl;
			if (param.first == "" ? true : ImGui::CollapsingHeader(param.first.c_str())) {
				ImGui::PushID(param.first.c_str());
				for (auto p : param.second) {
					//std::cout << p.first << "::"<<p.second << std::endl;
					if (uiFunctions.find(p.second->type) == uiFunctions.end()) {
						ImGui::Text("%s", p.first.c_str());
					}
					else {
						uiFunctions[p.second->type](*p.second);
					}
				}
				ImGui::PopID();
			}
		}
		ImGui::End();
	}
}

