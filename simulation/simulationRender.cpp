#include <glad/glad.h>
#include <GLFW/glfw3.h>
// include order for gl has to be left like this
#include "SPH.h"
#include "SPHRender.h"
#include <math.h>
#include "tinycolormap.h"
#include <chrono>
#include <iostream>
#include <sstream>
#include <vector>

#define WATCH(type, ns, id)\
  static type id = ParameterManager::instance().get<type>(std::string(#ns) +"." + #id);\
  if (ParameterManager::instance().get<type>(std::string(#ns) +"." + #id) != id) {\
      id = ParameterManager::instance().get<type>(std::string(#ns) +"." + #id);\
      dirty = true;\
  }


GLuint create1DTexture(v4* colorMap, int32_t elements) {
    GLuint textureId_;

    // generate the specified number of texture objects
    glGenTextures(1, &textureId_);
     assert(glGetError() == GL_NO_ERROR);

    // bind texture
    glBindTexture(GL_TEXTURE_1D, textureId_);
     assert(glGetError() == GL_NO_ERROR);

    // tells OpenGL how the data that is going to be uploaded is aligned
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
     assert(glGetError() == GL_NO_ERROR);

    glTexImage1D(
        GL_TEXTURE_1D, // Specifies the target texture. Must be GL_TEXTURE_1D or GL_PROXY_TEXTURE_1D.
        0, // Specifies the level-of-detail number. Level 0 is the base image level. Level n is the
           // nth mipmap reduction image.
        GL_RGBA32F, elements,
        0, // border: This value must be 0.
        GL_RGBA, GL_FLOAT, colorMap);
     assert(glGetError() == GL_NO_ERROR);

    // texture sampling/filtering operation.
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
     assert(glGetError() == GL_NO_ERROR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
     assert(glGetError() == GL_NO_ERROR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
     assert(glGetError() == GL_NO_ERROR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
     assert(glGetError() == GL_NO_ERROR);

    glBindTexture(GL_TEXTURE_1D, 0);
    // assert(glGetError() == GL_NO_ERROR);

    return textureId_;
}
#define STB_IMAGE_IMPLEMENTATION
#include <tools/stb_image.h>

GLuint createProgram(std::string vertexSource, std::string fragmentSource) {
    // Create an empty vertex shader handle
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);

    // Send the vertex shader source code to GL
    // Note that std::string's .c_str is NULL character terminated.
    const GLchar* source = (const GLchar*)vertexSource.c_str();
    glShaderSource(vertexShader, 1, &source, 0);

    // Compile the vertex shader
    glCompileShader(vertexShader);

    GLint isCompiled = 0;
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &isCompiled);
    if (isCompiled == GL_FALSE) {
        GLint maxLength = 0;
        glGetShaderiv(vertexShader, GL_INFO_LOG_LENGTH, &maxLength);

        // The maxLength includes the NULL character
        std::vector<GLchar> infoLog(maxLength);
        glGetShaderInfoLog(vertexShader, maxLength, &maxLength, &infoLog[0]);

        // We don't need the shader anymore.
        glDeleteShader(vertexShader);

        // Use the infoLog as you see fit.
        std::cerr << infoLog.data() << std::endl;

        // In this simple program, we'll just leave
        return -1;
    }

    // Create an empty fragment shader handle
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

    // Send the fragment shader source code to GL
    // Note that std::string's .c_str is NULL character terminated.
    source = (const GLchar*)fragmentSource.c_str();
    glShaderSource(fragmentShader, 1, &source, 0);

    // Compile the fragment shader
    glCompileShader(fragmentShader);

    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &isCompiled);
    if (isCompiled == GL_FALSE) {
        GLint maxLength = 0;
        glGetShaderiv(fragmentShader, GL_INFO_LOG_LENGTH, &maxLength);

        // The maxLength includes the NULL character
        std::vector<GLchar> infoLog(maxLength);
        glGetShaderInfoLog(fragmentShader, maxLength, &maxLength, &infoLog[0]);

        // We don't need the shader anymore.
        glDeleteShader(fragmentShader);
        // Either of them. Don't leak shaders.
        glDeleteShader(vertexShader);
        std::cerr << infoLog.data() << std::endl;

        // Use the infoLog as you see fit.

        // In this simple program, we'll just leave
        return -1;
    }

    // Vertex and fragment shaders are successfully compiled.
    // Now time to link them together into a program.
    // Get a program object.
    GLuint program = glCreateProgram();

    // Attach our shaders to our program
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);

    // Link our program
    glLinkProgram(program);

    // Note the different functions here: glGetProgram* instead of glGetShader*.
    GLint isLinked = 0;
    glGetProgramiv(program, GL_LINK_STATUS, (int*)&isLinked);
    if (isLinked == GL_FALSE) {
        GLint maxLength = 0;
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &maxLength);

        // The maxLength includes the NULL character
        std::vector<GLchar> infoLog(maxLength);
        glGetProgramInfoLog(program, maxLength, &maxLength, &infoLog[0]);
        std::cerr << infoLog.data() << std::endl;

        // We don't need the program anymore.
        glDeleteProgram(program);
        // Don't leak shaders either.
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);

        // Use the infoLog as you see fit.

        // In this simple program, we'll just leave
        return -1;
    }

    // Always detach shaders after a successful link.
    glDetachShader(program, vertexShader);
    glDetachShader(program, fragmentShader);
    return program;
}


void initRender() {

    glClearColor(1.f, 1.f, 1.f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT);
    glEnable(GL_POINT_SMOOTH);
    glPointSize(0.25f);
    glLineWidth(5.f);
    glMatrixMode(GL_PROJECTION);
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glClear(GL_COLOR_BUFFER_BIT);
}

#include <random>
#include "2DMath.h"


enum struct buffer_t {
    velocity, angularVelocity, acceleration, density, neighbors, UID, alpha, area, pressure1, pressure2, pressureBoundary, source, dpdt, rhoStar, predictedVelocity, pressureAcceleration, vorticity
};
enum struct renderMode_t {
    real, imag, realGradient, imagGradient, fractal
};


std::pair<float, float> updateField(float* data) {
    //static std::random_device rd;
    //static std::default_random_engine generator(rd());
    //static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    //std::cout << "Regenerating Data! " << std::endl;
    //for (int32_t i = 0; i < screenWidth * screenHeight; ++i) {
    //    data[i] = distribution(generator);
    //}

    auto modestr = ParameterManager::instance().get<std::string>("colorMap.renderMode");

    renderMode_t mode = renderMode_t::real;
    if (modestr == "real") mode = renderMode_t::real;
    if (modestr == "imag") mode = renderMode_t::imag;
    if (modestr == "real gradient") mode = renderMode_t::realGradient;
    if (modestr == "imag gradient") mode = renderMode_t::imagGradient;
    if (modestr == "fractal") mode = renderMode_t::fractal;
    {
    bool dirty = false;
    WATCH(int32_t, field, nx);
    WATCH(int32_t, field, ny);
    if(scalarFieldData == nullptr)
        scalarFieldData = (float*) malloc(sizeof(float) * nx * ny);
    if(dirty){
        scalarFieldData = (float*) realloc(scalarFieldData, sizeof(float) * nx * ny);
    }
    }
    auto hScale = ParameterManager::instance().get<scalar>("field.h");
    auto ny = ParameterManager::instance().get<int32_t>("field.ny");
    auto nx = ParameterManager::instance().get<int32_t>("field.nx");
    auto hi = (int32_t)::ceil(hScale);

    auto [domainMin, domainMax] = getDomain();

    #pragma omp parallel for
    for (int32_t y = 0; y < nx; ++y) {
        for (int32_t x = 0; x < ny; ++x) {
            scalar xPos = ((scalar) x) / ((scalar) nx - 1.) * (scalar) (domainMax.x() - domainMin.x()) + domainMin.x();
            scalar yPos = ((scalar) y) / ((scalar) ny - 1.) * (scalar) (domainMax.y() - domainMin.y()) + domainMin.y();

            auto val = evalFunction(cheb::complex(xPos, yPos));

            auto value = 0.;
            switch(mode){
case renderMode_t::real:value = val.real(); break;
case renderMode_t::imag:value = val.imag(); break;
case renderMode_t::realGradient:value = std::abs(val); break;
            }

            scalarFieldData[y * nx + x] = (float) value;
        }
    }

    auto min = *std::min_element(scalarFieldData, scalarFieldData + nx * ny);
    auto max = *std::max_element(scalarFieldData, scalarFieldData + nx * ny);
#pragma omp parallel for
    for (int32_t y = 0; y < screenHeight; ++y) {
        for (int32_t x = 0; x < screenWidth; ++x) {
            int32_t xi = std::clamp((int32_t)::floor(((scalar) x) / ((scalar) screenWidth) * (scalar) nx),0,nx - 1);
            int32_t yi = std::clamp((int32_t)::floor(((scalar) y) / ((scalar) screenHeight) * (scalar) ny),0,ny - 1);
            
            auto v = scalarFieldData[yi * nx + xi];
            v = (v - min) / (max - min);
            v = std::clamp(v, 0.f, 1.f);
            data[y * screenWidth + x] = v;
        }
    }
    return std::make_pair(min, max);

}

void renderField() {
    static GLuint textureID, vao, texID, program, minUni, maxUni, texUnit;
    static bool once = true;
    static float* data = new float[screenWidth * screenHeight];
    if(once){
        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        for (int32_t i = 0; i < screenWidth * screenHeight; ++i) {
            data[i] = 0.f;
        }
    auto vtxShader = R"(#version 150

// Input vertex data, different for all executions of this shader.
in vec3 vertexPosition_modelspace;

// Output data ; will be interpolated for each fragment.
out vec2 UV;

void main(){
	gl_Position =  vec4(vertexPosition_modelspace,1);
	UV = (vertexPosition_modelspace.xy+vec2(1,1))/2.0;
})";
    auto frgShader = R"(#version 150

in vec2 UV;

out vec3 color;

uniform sampler2D renderedTexture;
uniform sampler1D colorRamp;
uniform float min = 0.0;
uniform float max = 1.0;

void main(){
float v = texture( renderedTexture, UV ).r;
float rel = v;//(v - min) / (max - min);
float clamped = clamp(rel, 0.0, 1.0);
vec3 col = texture(colorRamp, clamped).xyz;
color = col;
	//color =  vec3(v,v,v);
	//color =  vec3(UV.xy,1.0);
})";
    GLuint quad_VertexArrayID;
    glGenVertexArrays(1, &quad_VertexArrayID);
    glBindVertexArray(quad_VertexArrayID);

    static const GLfloat g_quad_vertex_buffer_data[] = {
        -1.0f, -1.0f, 0.0f,
        1.0f, -1.0f, 0.0f,
        -1.0f,  1.0f, 0.0f,
        -1.0f,  1.0f, 0.0f,
        1.0f, -1.0f, 0.0f,
        1.0f,  1.0f, 0.0f,
    };
    // Create one OpenGL texture
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glGenTextures(1, &textureID);

    // "Bind" the newly created texture : all future texture functions will modify this texture
    glBindTexture(GL_TEXTURE_2D, textureID);

    // Give the image to OpenGL
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, screenWidth, screenHeight, 0, GL_RED, GL_FLOAT, data);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    GLuint quad_vertexbuffer;
    glGenBuffers(1, &quad_vertexbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, quad_vertexbuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_quad_vertex_buffer_data), g_quad_vertex_buffer_data, GL_STATIC_DRAW);

    // Create and compile our GLSL program from the shaders
    program = createProgram(vtxShader, frgShader);

    glUseProgram(program);
    texID = glGetUniformLocation(program, "renderedTexture");
    minUni = glGetUniformLocation(program, "min");
    maxUni = glGetUniformLocation(program, "max");
    //std::cout << "texID: " << texID << std::endl;
    //std::cout << "minUni: " << minUni << std::endl;
    //std::cout << "maxUni: " << maxUni << std::endl;

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, textureID);
    // Set our "myTextureSampler" sampler to use Texture Unit 0
    glUniform1i(texID, 1);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, quad_vertexbuffer);
    glVertexAttribPointer(
        0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
        3,                  // size
        GL_FLOAT,           // type
        GL_FALSE,           // normalized?
        0,                  // stride
        (void*)0            // array buffer offset
    );
    int32_t image_width = 1024;
    int32_t image_height = 1024;
    v4* img = new v4[1024];
    v4* color_map = nullptr;
    for (int32_t it = 0; it < 1024; ++it)
        img[it] = v4{ (float)it / (float)1024 * 255.f,(float)it / (float)1024 * 255.f,(float)it / (float)1024 * 255.f, 255.f };
    std::string file_name = "viridis.png";
    if (std::filesystem::exists(file_name)) {
        std::cout << "Loading " << file_name << std::endl;
        unsigned char* image_data = stbi_load(file_name.c_str(), &image_width, &image_height, NULL, 4);
        delete[] img;
        img = new v4[image_width];
        for (int32_t it = 0; it < image_width; ++it) {
            img[it] = v4{
            (float)image_data[it * 4 + 0],
            (float)image_data[it * 4 + 1],
            (float)image_data[it * 4 + 2], 255.f };
            //std::cout << it << ": [ " << img[it].x << " " << img[it].y << " " << img[it].z <<
            //    " ]" << std::endl;
        }
        //img = QImage(QString::fromStdString(file_name));
        //img.load(QString(file_name.c_str()));
        //std::cout << image_width << " : " << image_height << std::endl;
    }
    //catch (...) {}
    color_map = (v4*)realloc(color_map, sizeof(v4) * (image_width));
    for (int32_t it = 0; it < image_width; ++it) {
        color_map[it] = v4{ (float)(img[it].x) / 256.f, (float)(img[it].y) / 256.f,
                                  (float)(img[it].z) / 256.f, 1.f };
        //std::cout << color_map[it].x() << " : " << color_map[it].y() << " : " << color_map[it].z() << std::endl;
          //if(it == img.width()-1)
          //	color_map[it + 1] = QVector4D{ (float)(col.red()) / 256.f, (float)(col.green()) / 256.f,
          //	(float)(col.blue()) / 256.f, 1.f };
    }
    auto color_map_elements = image_width;

    //std::cout << "color_map_elements " << color_map_elements << std::endl;
    texUnit = create1DTexture(color_map, color_map_elements);
    GLuint samplerLocation = glGetUniformLocation(program, "colorRamp");
    //std::cout << "colorRamp: " << samplerLocation << std::endl;
    glUseProgram(program);
    glUniform1i(samplerLocation, 0);
    glActiveTexture(GL_TEXTURE0 + 0);
    glBindTexture(GL_TEXTURE_1D, texUnit);
    //std::cout << " unit " << texUnit << " ->  " << samplerLocation << " : " << 0 << std::endl;

    glUseProgram(0);
    glBindVertexArray(0);
}

static std::string colorMap = "";
//std::cout << colorMap << " : " << ParameterManager::instance().get<std::string>("colorMap.colorMap") << std::endl; 
if(ParameterManager::instance().get<std::string>("colorMap.colorMap") != colorMap){
    colorMap = ParameterManager::instance().get<std::string>("colorMap.colorMap");
    int32_t image_width = 1024;
    int32_t image_height = 1024;
    v4* img = new v4[1024];
    v4* color_map = nullptr;
    for (int32_t it = 0; it < 1024; ++it)
        img[it] = v4{ (float)it / (float)1024 * 255.f,(float)it / (float)1024 * 255.f,(float)it / (float)1024 * 255.f, 255.f };
    // std::string file_name = std::string("cfg/") + colorMap + ".png";
    
    std::string file_name = resolveFile(std::string("cfg/") + ParameterManager::instance().get<std::string>("colorMap.colorMap") + ".png").string();

    if (std::filesystem::exists(file_name)) {
        std::cout << "Loading " << file_name << std::endl;
        unsigned char* image_data = stbi_load(file_name.c_str(), &image_width, &image_height, NULL, 4);
        delete[] img;
        img = new v4[image_width];
        for (int32_t it = 0; it < image_width; ++it) {
            img[it] = v4{
            (float)image_data[it * 4 + 0],
            (float)image_data[it * 4 + 1],
            (float)image_data[it * 4 + 2], 255.f };
            //std::cout << it << ": [ " << img[it].x << " " << img[it].y << " " << img[it].z <<
            //    " ]" << std::endl;
        }
        //img = QImage(QString::fromStdString(file_name));
        //img.load(QString(file_name.c_str()));
        //std::cout << image_width << " : " << image_height << std::endl;
    }
    //catch (...) {}
    color_map = (v4*)realloc(color_map, sizeof(v4) * (image_width));
    for (int32_t it = 0; it < image_width; ++it) {
        color_map[it] = v4{ (float)(img[it].x) / 256.f, (float)(img[it].y) / 256.f,
                                  (float)(img[it].z) / 256.f, 1.f };
        //std::cout << color_map[it].x() << " : " << color_map[it].y() << " : " << color_map[it].z() << std::endl;
          //if(it == img.width()-1)
          //	color_map[it + 1] = QVector4D{ (float)(col.red()) / 256.f, (float)(col.green()) / 256.f,
          //	(float)(col.blue()) / 256.f, 1.f };
    }
    auto color_map_elements = image_width;

    //std::cout << "color_map_elements " << color_map_elements << std::endl;
    texUnit = create1DTexture(color_map, color_map_elements);
}


    bool dirty = once;

  WATCH(std::string, colorMap, renderMode);
  WATCH(int32_t, field, nx);
  WATCH(int32_t, field, ny);
  WATCH(scalar, domain, xmin);
  WATCH(scalar, domain, xmax);
  WATCH(scalar, domain, ymin);
  WATCH(scalar, domain, ymax);

    if (dirty) {
        auto [min,max] = updateField(data);
        //std::cout << buff << " -> " << vMode << std::endl;
        //std::cout << "New min and max: " << min << " x " << max << std::endl;
        ParameterManager::instance().get<scalar>("field.min") = min;
        ParameterManager::instance().get<scalar>("field.max") = max;
        glUseProgram(program);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, screenWidth, screenHeight, GL_RED, GL_FLOAT, data);
        glUseProgram(0);
    }
    glUseProgram(program);
    //std::cout << textureID << " -> " << texUnit << std::endl;
    glBindVertexArray(vao);
    glUniform1f(minUni, ParameterManager::instance().get<scalar>("field.min"));
    glUniform1f(maxUni, ParameterManager::instance().get<scalar>("field.max"));
    //GLuint samplerLocation = glGetUniformLocation(program, "colorRamp");
    //glUniform1i(samplerLocation, 0);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, textureID);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_1D, texUnit);
    glDrawArrays(GL_TRIANGLES, 0, 6); // 2*3 indices starting at 0 -> 2 triangles
    glBindVertexArray(0);
    glUseProgram(0);
    once = false;
}

void render() {
    glPointSize(0.001f);
    glPointSize(3.f);
    // used to find the proper scaling for color mapping
    //auto velmap = [](const Particle p) -> scalar { return p.vel.norm(); };
    //auto rhomap = [](const Particle p) -> scalar { return p.rho; };
    //auto neighmap = [](const Particle p) -> scalar {
    //    return (scalar) p.neighbors.size(); 
    //};

    // reset the openGL state
    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glOrtho(0, 1.0, 0, 1.0, 0, 1);
    //glOrtho(40, 60, 3, 13, 0, 1);

        renderField();

    glColor4f(0.8f, 0.f, 0.f, 1);
}