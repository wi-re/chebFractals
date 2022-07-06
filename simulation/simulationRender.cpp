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
      std::cout << "Parameter " << std::string(#ns) << "." << std::string(#id) << " changed from " << id << " to " << ParameterManager::instance().get<type>(std::string(#ns) +"." + #id) << std::endl;\
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
    real, imag, realGradient, imagGradient, fractal, abs, absGradient, f2, Jx, Jy
};

#include <mutex>
enum struct optimizationMethod{
    newton, halley, gradientDescent, newtonOptimizer, newtonOptimizerHessian, Adam,gradientDescentAdaptive,adaGrad,BFGS
};
auto getMethod(std::string stringMethod){    
    if(stringMethod == "newton")
        return optimizationMethod::newton;
    if(stringMethod == "halley")
        return optimizationMethod::halley;
    if(stringMethod == "gradientDescent")
        return optimizationMethod::gradientDescent;
    if(stringMethod == "gradientDescentAdaptive")
        return optimizationMethod::gradientDescentAdaptive;
    if(stringMethod == "newton optimizer")
        return optimizationMethod::newtonOptimizer;
    if(stringMethod == "newton hessian")
        return optimizationMethod::newtonOptimizerHessian;
    if(stringMethod == "adam")
        return optimizationMethod::Adam;
        if(stringMethod=="adaGrad")
        return optimizationMethod::adaGrad;
        if(stringMethod=="BFGS")
        return optimizationMethod::BFGS;
    return optimizationMethod::newton;
}
enum optimizationState : int32_t{
    converged, diverged, unstable
};

auto evalSQFunction(Eigen::Vector2d location){
    auto [f, J, H] = evalSquarePolynomial(cheb::complex(location.x(), location.y()));
    return f;
}
using v2 = Eigen::Vector2d;
using mat2 = Eigen::Matrix2d;
auto finiteGrad(Eigen::Vector2d p){
    auto [f, J, H] = evalSquarePolynomial(cheb::complex(p.x(), p.y()));
    return v2(J.dfdx, J.dfdy);


    auto h = std::cbrt(DBL_EPSILON);

    auto dx = (evalSQFunction(v2{p.x() + h,p.y()})) - evalSQFunction(v2{p.x() - h,p.y()}) / (2. * h);
    auto dy = (evalSQFunction(v2{p.x(),p.y() + h})) - evalSQFunction(v2{p.x(),p.y() - h}) / (2. * h);

    return v2(dx,dy);
}

auto lineSearchInv(Eigen::Vector2d x, Eigen::Vector2d p){
    auto f = [](Eigen::Vector2d x){
        return evalSQFunction(x);
    };
    auto nabla = [](Eigen::Vector2d x){
        return finiteGrad(x);
    };

    auto alpha = ParameterManager::instance().get<scalar>("BFGS.alpha");
    static auto& c1 = ParameterManager::instance().get<scalar>("BFGS.c1");
    static auto& c2 = ParameterManager::instance().get<scalar>("BFGS.c2");
    static auto& epsp = ParameterManager::instance().get<scalar>("BFGS.eps");
    auto eps = std::pow(10., epsp);
    // auto alpha = -1.;
    // auto c1 = 1e-4;
    // auto c2 = 0.9;
    int32_t i = 0;

    bool armijo = f(x + alpha * p) <= f(x) + c1 * alpha * p.dot(nabla(x));
    bool curvature = -p.dot(nabla(x + alpha * p)) <= - c2 * p.dot(nabla(x));

    while((!armijo || !curvature) && std::abs(alpha) >= eps){
        alpha *= 0.5;
        auto x_new = x + alpha * p;
        auto nabla_new = finiteGrad(x_new);
        armijo = f(x + alpha * p) <= f(x) + c1 * alpha * p.dot(nabla(x));
        curvature = -p.dot(nabla(x + alpha * p)) <= - c2 * p.dot(nabla(x));
        // printf(
        //     "\t[%03d]: [%g %g] + %g * [%g %g] = [%g %g] -> [%g %g] -> %g >= %g || %g <= %g\n",
        //     ++i, x.x(), x.y(), alpha, p.x(), p.y(), x_new.x(), x_new.y(), nabla_new.x(), nabla_new.y(),
        //     f(x + alpha * p) , f(x) + c1 * alpha * p.dot(nabla(x)),
        //     -p.dot(nabla(x + alpha * p)) , - c2 * p.dot(nabla(x)));

    }
    return alpha;
}
auto lineSearch(Eigen::Vector2d x, Eigen::Vector2d p){
    auto f = [](Eigen::Vector2d x){
        return evalSQFunction(x);
    };
    auto nabla = [](Eigen::Vector2d x){
        return finiteGrad(x);
    };
    auto alpha = ParameterManager::instance().get<scalar>("BFGS.alpha");
    static auto& c1 = ParameterManager::instance().get<scalar>("BFGS.c1");
    static auto& c2 = ParameterManager::instance().get<scalar>("BFGS.c2");
    static auto& epsp = ParameterManager::instance().get<scalar>("BFGS.eps");
    auto eps = std::pow(10., epsp);

    int32_t i = 0;

    bool armijo = f(x + alpha * p) <= f(x) + c1 * alpha * p.dot(nabla(x));
    bool curvature = -p.dot(nabla(x + alpha * p)) <= - c2 * p.dot(nabla(x));

    while((!armijo || !curvature) && alpha >= eps){
        alpha *= 0.5;
        auto x_new = x + alpha * p;
        auto nabla_new = finiteGrad(x_new);
        armijo = f(x + alpha * p) <= f(x) + c1 * alpha * p.dot(nabla(x));
        curvature = -p.dot(nabla(x + alpha * p)) <= - c2 * p.dot(nabla(x));
        // printf(
        //     "\t[%03d]: [%g %g] + %g * [%g %g] = [%g %g] -> [%g %g] -> %g >= %g || %g <= %g\n",
        //     ++i, x.x(), x.y(), alpha, p.x(), p.y(), x_new.x(), x_new.y(), nabla_new.x(), nabla_new.y(),
        //     f(x + alpha * p) , f(x) + c1 * alpha * p.dot(nabla(x)),
        //     -p.dot(nabla(x + alpha * p)) , - c2 * p.dot(nabla(x)));

    }
    if(alpha < eps)
        return lineSearchInv(x,p);
    return alpha;
}

auto BFGS(cheb::complex location){
    std::vector<cheb::complex> positions, values, steps;
    static auto& stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
    static auto& learningRatef = ParameterManager::instance().get<scalar>("field.learningRate");

    auto x0 = v2{location.real(), location.imag()};
    
    auto [f, J, H] = evalSquarePolynomial(cheb::complex(location.real(), location.imag()));

    mat2 Bk;
    // Bk << H.d2fdx2, H.d2fdxy, H.d2fdyx, H.d2fdy2;
    Bk << 1.,0.,0.,1.;
    Bk = Bk.inverse();
    v2 x = x0;
    int32_t it = 2;

    positions.push_back(cheb::complex(x.x(), x.y()));
    values.push_back(evalSQFunction(x));
    auto step = cheb::complex(0.,0.);
    steps.push_back(step);


    for(int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i){    
        it += 1;
        v2 p = -Bk * finiteGrad(x);
        //printf("[%03d]: %g + %gi => Starting Linesearch\n",it, p.x(), p.y());

        auto a = lineSearch(x,p);
        auto s = a * p ;
        //printf("[%03d]: a = %g : [%g + %gi]\n",it, a, s.x(), s.y());
        v2 x_new = x + s;

        auto y = finiteGrad(x_new) - finiteGrad(x);

        scalar sTy = s.dot(y);
        Eigen::RowVector2d yTB = y.transpose() * Bk;
        scalar yTBy = yTB * y;
        mat2 ssT = s * s.transpose();   
        mat2 leftTerm = (sTy + yTBy) / (sTy * sTy) * ssT;

        mat2 BysT = (Bk * y) * s.transpose();
        mat2 syTB = (s * y.transpose()) * Bk;
        mat2 rightTerm = (BysT + syTB) / (sTy);

        //mat2 Bkp = Bk + leftTerm - rightTerm;

        mat2 Id = mat2::Identity();
        scalar rhok = 1. / y.dot(s);
        mat2 left = Id - rhok * s * y.transpose();
        mat2 right = Id - rhok * y * s.transpose();
        mat2 righter = rhok * s * s.transpose();
        mat2 Bkp = left * Bk * right + righter;




        //printf("[%03d]: f(%g + %gi) = %g : [%g  + %gi] -> [%g %g] [%g %g]\n",it, x_new.x(), x_new.y(),evalSQFunction(x_new), finiteGrad(x_new).x(), finiteGrad(x_new).y(), Bkp(0,0), Bkp(0,1), Bkp(1,0), Bkp(1,1));

        // mat2 Bkp = Bk + (s.transpose() * y + y.transpose() * Bk * y) *(s * s.transpose()) / (s.transpose() * y).squaredNorm() - (Bk * y * s.transpose() + s * y.transpose() * Bk) / (s.transpose() * y)[0];

        
        positions.push_back(cheb::complex(x.x(), x.y()));
        values.push_back(evalSQFunction(x));
        auto step = cheb::complex(s.x(), s.y());
        steps.push_back(step);

        if(std::abs(step) < 1e-14)
            return std::make_tuple(optimizationState::converged, location, positions, values, steps);
        if(location != location || step != step)
            return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
        location += step;
        x = x_new;

        Bk = Bkp;
    }
    return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}




auto newtonsMethod(cheb::complex location){
    std::vector<cheb::complex> positions, values, steps;
    static auto& stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
    static auto& learningRatef = ParameterManager::instance().get<scalar>("field.learningRate");
    auto learningRate = std::pow(10.,learningRatef);

    for(int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i){      
        auto fx = evalFunction(location);
        auto dx = evalDerivative(location);
        auto step = -fx / dx;
        positions.push_back(location);
        values.push_back(fx);
        steps.push_back(step);

        if(std::abs(step) < 1e-12)
            return std::make_tuple(optimizationState::converged, location, positions, values, steps);
        if(location != location || step != step)
            return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
        location += step;
    }
    return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}
auto newtonsMethodOptimizer(cheb::complex location){
    std::vector<cheb::complex> positions, values, steps;
    static auto& stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
    static auto& learningRatef = ParameterManager::instance().get<scalar>("field.learningRate");
    auto learningRate = std::pow(10.,learningRatef);

    for(int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i){      
        auto [f, J, H] = evalSquarePolynomial(location);

        std::complex a = H.d2fdx2;
        std::complex b = H.d2fdxy;
        std::complex c = H.d2fdyx;
        std::complex d = H.d2fdy2;
        std::complex T = a + d;
        std::complex D = a*d-b*c;
        std::complex l1 = T/2. + std::sqrt(T * T / 4. - D);      
        std::complex l2 = T/2. - std::sqrt(T * T / 4. - D);      
        auto det = H.d2fdx2 * H.d2fdy2 - H.d2fdxy * H.d2fdyx;

        Hessian H_1{ H.d2fdy2 / det, -H.d2fdxy / det, -H.d2fdyx / det, H.d2fdx2/det};

        auto prodx = H_1.d2fdx2 * J.dfdx + H_1.d2fdxy * J.dfdy;
        auto prody = H_1.d2fdyx * J.dfdx + H_1.d2fdy2 * J.dfdy;

        auto dx = cheb::complex(prodx, prody);
        auto step = -dx * 1.;

        positions.push_back(location);
        values.push_back(f);
        steps.push_back(step);        if(std::abs(step) < 1e-12)
            return std::make_tuple(optimizationState::converged, location, positions, values, steps);
        if(location != location || step != step)
            return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
        location += step;
    }
    return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}
auto newtonsMethodOptimizerHessian(cheb::complex location){
    std::vector<cheb::complex> positions, values, steps;
    static auto& stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
    static auto& learningRatef = ParameterManager::instance().get<scalar>("field.learningRate");
    auto learningRate = std::pow(10.,learningRatef);

    for(int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i){      
        auto [f, J, H] = evalSquarePolynomial(location);
        
        std::complex a = H.d2fdx2;
        std::complex b = H.d2fdxy;
        std::complex c = H.d2fdyx;
        std::complex d = H.d2fdy2;
        std::complex T = a + d;
        std::complex D = a*d-b*c;
        std::complex l1 = T/2. + std::sqrt(T * T / 4. - D);      
        std::complex l2 = T/2. - std::sqrt(T * T / 4. - D);         
        auto l = std::min(l1.real(), l2.real());
        if(l < 0.){
            H.d2fdx2 -= l * 1.5;
            H.d2fdy2 -= l * 1.5;
        }

        auto det = H.d2fdx2 * H.d2fdy2 - H.d2fdxy * H.d2fdyx;

        Hessian H_1{ H.d2fdy2 / det, -H.d2fdxy / det, -H.d2fdyx / det, H.d2fdx2/det};

        auto prodx = H_1.d2fdx2 * J.dfdx + H_1.d2fdxy * J.dfdy;
        auto prody = H_1.d2fdyx * J.dfdx + H_1.d2fdy2 * J.dfdy;

        auto dx = cheb::complex(prodx, prody);
        auto step = -dx * 1.;

        positions.push_back(location);
        values.push_back(f);
        steps.push_back(step);
                if(std::abs(step) < 1e-12)
            return std::make_tuple(optimizationState::converged, location, positions, values, steps);
        if(location != location || step != step)
            return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
        location += step;
    }
    return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}
auto gradientDescentAdaptive(cheb::complex location){
    std::vector<cheb::complex> positions, values, steps;
    static auto& stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
    static auto& learningRatef = ParameterManager::instance().get<scalar>("field.learningRate");
    auto learningRate = std::pow(10.,learningRatef);

    auto [f, J, H] = evalSquarePolynomial(location);
    auto dxPrior = cheb::complex(J.dfdx,J.dfdy);
    auto locationPrior = location;

    auto step = -dxPrior * learningRate;
    positions.push_back(location);
    values.push_back(f);
    steps.push_back(step);


    location += step;

    for(int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i){   
        auto [f, J, H] = evalSquarePolynomial(location);

        auto dx = cheb::complex(J.dfdx,J.dfdy);
        // dx = cheb::complex(prodx, prody);

        auto diff = location - locationPrior;
        auto diffGrad = dx - dxPrior;

        auto gamma =  (diff.real()     * diffGrad.real() + diff.imag()     * diffGrad.imag()) / 
                      (diffGrad.real() * diffGrad.real() + diffGrad.imag() * diffGrad.imag() +  1e-12);



        auto step = -gamma * dx;


    
        positions.push_back(location);
        values.push_back(f);
        steps.push_back(step);
               if(std::abs(step) < 1e-12)
            return std::make_tuple(optimizationState::converged, location, positions, values, steps);
        if(location != location || step != step)
            return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
        locationPrior = location;
        location += step;
        dxPrior = dx;
    }
    return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}

auto adaGrad(cheb::complex location){
    std::vector<cheb::complex> positions, values, steps;
    static auto& stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
    static auto& learningRatef = ParameterManager::instance().get<scalar>("field.learningRate");
    auto learningRate = std::pow(10.,learningRatef);

    auto sumx = 0.;
    auto sumy = 0.;

    for(int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i){   
        auto [f, J, H] = evalSquarePolynomial(location);

        auto dx = cheb::complex(J.dfdx,J.dfdy);
        // dx = cheb::complex(prodx, prody);

        sumx += J.dfdx * J.dfdx;
        sumy += J.dfdy * J.dfdy;

        auto alphax = learningRate / (1e-8 + std::sqrt(sumx));
        auto alphay = learningRate / (1e-8 + std::sqrt(sumy));
        auto step = - cheb::complex(J.dfdx * alphax, J.dfdy * alphay);


        // auto step = -dx * learningRate;


    
        positions.push_back(location);
        values.push_back(f);
        steps.push_back(step);
               if(std::abs(step) < 1e-12)
            return std::make_tuple(optimizationState::converged, location, positions, values, steps);
        if(location != location || step != step)
            return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
        location += step;
    }
    return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}
auto gradientDescent(cheb::complex location){
    std::vector<cheb::complex> positions, values, steps;
    static auto& stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
    static auto& learningRatef = ParameterManager::instance().get<scalar>("field.learningRate");
    auto learningRate = std::pow(10.,learningRatef);

    for(int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i){   
        auto [f, J, H] = evalSquarePolynomial(location);

        auto dx = cheb::complex(J.dfdx,J.dfdy);
        // dx = cheb::complex(prodx, prody);

        auto step = -dx * learningRate;


    
        positions.push_back(location);
        values.push_back(f);
        steps.push_back(step);
               if(std::abs(step) < 1e-12)
            return std::make_tuple(optimizationState::converged, location, positions, values, steps);
        if(location != location || step != step)
            return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
        location += step;
    }
    return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}
auto halleysMethod(cheb::complex location){
    std::vector<cheb::complex> positions, values, steps;
    static auto& stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
    static auto& learningRatef = ParameterManager::instance().get<scalar>("field.learningRate");
    auto learningRate = std::pow(10.,learningRatef);

    for(int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i){        
        auto [f, J, H] = evalPolynomial(location);
        auto step = -(2. * f * J) / (2. * J * J - f *H); 


        positions.push_back(location);
        values.push_back(f);
        steps.push_back(step);

               if(std::abs(step) < 1e-12)
            return std::make_tuple(optimizationState::converged, location, positions, values, steps);
        if(location != location || step != step)
            return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
        location += step;
    }
    return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}

auto adam(cheb::complex location){
    std::vector<cheb::complex> positions, values, steps;
    static auto& stepLimit = ParameterManager::instance().get<int32_t>("field.steps");
    static auto& alphaP = ParameterManager::instance().get<scalar>("adam.alpha");
    static auto& beta1 = ParameterManager::instance().get<scalar>("adam.beta1");
    static auto& beta2 = ParameterManager::instance().get<scalar>("adam.beta2");
    static auto& epsP = ParameterManager::instance().get<scalar>("adam.eps");

    auto alpha = std::pow(10., alphaP);
    auto eps = std::pow(10., epsP);

    auto t = 0.0;
    auto m_x = 0.;
    auto v_x = 0.;

    auto m_y = 0.;
    auto v_y = 0.;

    for(int32_t i = 0; stepLimit == 0 ? true : i < stepLimit; ++i){   
        auto [f, J, H] = evalSquarePolynomial(location);

        t = t + 1.;
        // get gradients
        auto g_t = J;



        m_x = beta1 * m_x + (1. - beta1) * g_t.dfdx;
        v_x = beta2 * v_x + (1. - beta2) * g_t.dfdx * g_t.dfdx; 
        auto m_x_corrected = m_x / (1. - std::pow(beta1, t));
        auto v_x_corrected = v_x / (1. - std::pow(beta2, t));
        auto step_x = - alpha * m_x_corrected / (std::sqrt(v_x_corrected) + eps);


        m_y = beta1 * m_y + (1. - beta1) * g_t.dfdy;
        v_y = beta2 * v_y + (1. - beta2) * g_t.dfdy * g_t.dfdy; 
        auto m_y_corrected = m_y / (1. - std::pow(beta1, t));
        auto v_y_corrected = v_y / (1. - std::pow(beta2, t));
        auto step_y = - alpha * m_y_corrected / (std::sqrt(v_y_corrected) + eps);


        auto step = cheb::complex(step_x, step_y);


    
        positions.push_back(location);
        values.push_back(f);
        steps.push_back(step);
               if(std::abs(step) < 1e-12)
            return std::make_tuple(optimizationState::converged, location, positions, values, steps);
        if(location != location || step != step)
            return std::make_tuple(optimizationState::diverged, location, positions, values, steps);
        location += step;
    }
    return std::make_tuple(optimizationState::unstable, location, positions, values, steps);
}
auto optimize(cheb::complex location, optimizationMethod method){
    switch(method){
        case optimizationMethod::newton: return newtonsMethod(location);
        case optimizationMethod::newtonOptimizer: return newtonsMethodOptimizer(location);
        case optimizationMethod::newtonOptimizerHessian: return newtonsMethodOptimizerHessian(location);
        case optimizationMethod::gradientDescent: return gradientDescent(location);
        case optimizationMethod::gradientDescentAdaptive: return gradientDescentAdaptive(location);
        case optimizationMethod::halley: return halleysMethod(location);
        case optimizationMethod::adaGrad: return adaGrad(location);
        case optimizationMethod::Adam: return adam(location);
        case optimizationMethod::BFGS: return BFGS(location);
        default: return newtonsMethod(location);
    }
}


struct progressBarLocal {
    std::chrono::high_resolution_clock::time_point start, lastTime;
    int32_t frames, lastFrame = 0;

    progressBarLocal(int32_t targetFrames) :frames(targetFrames) {
        start = std::chrono::high_resolution_clock::now();
    };

    void update(int32_t current, int32_t updateRate) {
        std::ios cout_state(nullptr);
        cout_state.copyfmt(std::cout);
        if (current >= lastFrame + updateRate) {
            auto progress = ((double)current) / ((double)frames);
            auto now = std::chrono::high_resolution_clock::now();
            lastTime = now;
            int barWidth = 70;
            std::cout << "Rendering " << std::setw(4) << current;
            if (frames != -1)
                std::cout << "/" << std::setw(4) << frames;
            std::cout << " [";
            int pos = barWidth * progress;
            for (int i = 0; i < barWidth; ++i) {
                if (i < pos)
                    std::cout << "=";
                else if (i == pos)
                    std::cout << ">";
                else
                    std::cout << " ";
            }
            std::cout << "] " << std::setw(3) << int(progress * 100.0) << " ";
            auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
            if (dur.count() < 100 || progress < 1e-3f) {
                std::cout << " ---/---s  ";
            }
            else {
                auto totalTime =
                    ((float)std::chrono::duration_cast<std::chrono::microseconds>(now - start).count()) / 1000.f / 1000.f;
                std::cout << std::fixed << std::setprecision(0) << " " << std::setw(3) << totalTime << "/" << std::setw(3)
                    << (totalTime / progress) << "s  ";
            }


        }
        std::cout << "\r";
        std::cout.flush();
        std::cout.copyfmt(cout_state);
    }
    void update(int32_t x, int32_t y, int32_t updateRate) {
        std::ios cout_state(nullptr);
        cout_state.copyfmt(std::cout);
        int32_t current = x + y * 1920;
        if (current >= lastFrame + updateRate) {
            auto progress = ((double)current) / ((double)frames);
            auto now = std::chrono::high_resolution_clock::now();
            lastTime = now;
            int barWidth = 70;
            std::cout << "Rendering " << std::setw(5) << x << std::setw(5) << y;
            if (frames != -1)
                std::cout << "/" << std::setw(5) << 1920 << std::setw(5) << 1080;
            std::cout << " [";
            int pos = barWidth * progress;
            for (int i = 0; i < barWidth; ++i) {
                if (i < pos)
                    std::cout << "=";
                else if (i == pos)
                    std::cout << ">";
                else
                    std::cout << " ";
            }
            std::cout << "] " << std::setw(3) << int(progress * 100.0) << " ";
            auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
            if (dur.count() < 100 || progress < 1e-3f) {
                std::cout << " ---/---s  ";
            }
            else {
                auto totalTime =
                    ((float)std::chrono::duration_cast<std::chrono::microseconds>(now - start).count()) / 1000.f / 1000.f;
                std::cout << std::fixed << std::setprecision(0) << " " << std::setw(3) << totalTime << "/" << std::setw(3)
                    << (totalTime / progress) << "s  ";
            }


        }
        std::cout << "\r";
        std::cout.flush();
        std::cout.copyfmt(cout_state);
    }
    void end() {
        std::cout << std::endl;
    }
};

std::pair<float, float> updateField(float* data, float* angular, float* radial) {
    //static std::random_device rd;
    //static std::default_random_engine generator(rd());
    //static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    //std::cout << "Regenerating Data! " << std::endl;
    //for (int32_t i = 0; i < screenWidth * screenHeight; ++i) {
    //    data[i] = distribution(generator);
    //}

auto learningRate = std::pow(10., ParameterManager::instance().get<scalar>("field.learningRate"));
auto stringMethod = ParameterManager::instance().get<std::string>("field.method");
optimizationMethod method = getMethod(stringMethod);

    static scalar prevClusterEps = ParameterManager::instance().get<scalar>("field.clusterEpsilon");
    static auto& curClusterEps = ParameterManager::instance().get<scalar>("field.clusterEpsilon");
    scalar clusterEpsilon = std::pow(10., curClusterEps);
    static std::vector<cheb::complex> clusterCenters;
    if(prevClusterEps!=curClusterEps){
        clusterCenters.clear();
        prevClusterEps = curClusterEps;
    }
    static scalar prevThreshold = ParameterManager::instance().get<scalar>("field.threshold");
    static auto& curThreshold = ParameterManager::instance().get<scalar>("field.threshold");
    scalar clusterThreshold = std::pow(10., curThreshold);
    if(prevThreshold!=curThreshold){
        prevThreshold = curThreshold;
    }

    auto modestr = ParameterManager::instance().get<std::string>("colorMap.renderMode");
    static cheb::complex* vectorFieldData;
    static cheb::complex* functionFieldData;
    static cheb::complex* inputFieldData;
    static int32_t* iterationFieldData;
    renderMode_t mode = renderMode_t::real;
    if (modestr == "real") mode = renderMode_t::real;
    if (modestr == "imag") mode = renderMode_t::imag;
    if (modestr == "abs") mode = renderMode_t::abs;
    if (modestr == "real gradient") mode = renderMode_t::realGradient;
    if (modestr == "imag gradient") mode = renderMode_t::imagGradient;
    if (modestr == "abs gradient") mode = renderMode_t::absGradient;
    if (modestr == "fractal") mode = renderMode_t::fractal;
    if (modestr == "f2") mode = renderMode_t::f2;
    if (modestr == "Jx") mode = renderMode_t::Jx;
    if (modestr == "Jy") mode = renderMode_t::Jy;
    {
    bool dirty = false;
    WATCH(int32_t, field, nx);
    WATCH(int32_t, field, ny);
    if(scalarFieldData == nullptr)
        scalarFieldData = (float*) malloc(sizeof(float) * nx * ny);
    if(iterationFieldData == nullptr)
        iterationFieldData = (int32_t*) malloc(sizeof(int32_t) * nx * ny);
    if(vectorFieldData == nullptr)
        vectorFieldData = (cheb::complex*) malloc(sizeof(vec) * nx * ny);
    if(inputFieldData == nullptr)
        inputFieldData = (cheb::complex*) malloc(sizeof(vec) * nx * ny);
    if(functionFieldData == nullptr)
        functionFieldData = (cheb::complex*) malloc(sizeof(vec) * nx * ny);

        dataWidth = nx;
        dataHeight = ny;
    if(dirty){
        int32_t ny = ParameterManager::instance().get<int32_t>("field.ny");
        int32_t nx = ParameterManager::instance().get<int32_t>("field.nx");
        auto total = sizeof(float) * nx * ny;
        std::cout << "Allocating " << sizeof(float) << " x " << nx << " x " << ny << "B -> " << total << std::endl;
        scalarFieldData = (float*) realloc(scalarFieldData, total);
        iterationFieldData = (int32_t*) realloc(scalarFieldData, total);
        vectorFieldData = (cheb::complex*) realloc(vectorFieldData, total * 4);
        inputFieldData = (cheb::complex*) realloc(inputFieldData, total * 4);
        functionFieldData = (cheb::complex*) realloc(functionFieldData, total * 4);
        dataWidth = nx;
        dataHeight = ny;
    }
    }
    auto hScale = ParameterManager::instance().get<scalar>("field.h");
    auto ny = ParameterManager::instance().get<int32_t>("field.ny");
    auto nx = ParameterManager::instance().get<int32_t>("field.nx");
    auto hi = (int32_t)::ceil(hScale);
    bool clustering = ParameterManager::instance().get<bool>("field.clustering");
    bool cycles = ParameterManager::instance().get<bool>("field.cycles");

    auto [domainMin, domainMax] = getDomain();

    // for(int32_t it = 0; it < nx * ny; ++it){{
    //     int32_t y = it % nx;
    //     int32_t x = it / nx;
    // }
    std::cout << "Starting Visualization" << std::endl;
    progressBarLocal pbl{nx*ny};
    std::atomic<int32_t> iterationCounter = 0;

    #pragma omp parallel for schedule(dynamic, 1)
    for(int32_t it = 0; it < nx * ny; ++it){{
        int32_t x = it % nx;
        int32_t y = it / nx;
            scalar xPos = ((scalar) x) / ((scalar) nx - 1.) * (scalar) (domainMax.x() - domainMin.x()) + domainMin.x();
            scalar yPos = ((scalar) y) / ((scalar) ny - 1.) * (scalar) (domainMax.y() - domainMin.y()) + domainMin.y();



            if(mode != renderMode_t::fractal){
                    auto fx = evalFunction(cheb::complex(xPos, yPos));
                    auto dx = evalDerivative(cheb::complex(xPos, yPos));
                    auto [f,J,H] = evalSquarePolynomial(cheb::complex(xPos, yPos));
                    //auto d2x = evalFunction(cheb::complex(xPos, yPos));

                    auto value = 0.;
                    switch(mode){
                        case renderMode_t::real:value = fx.real(); break;
                        case renderMode_t::imag:value = fx.imag(); break;
                        case renderMode_t::abs: value = std::abs(fx); break;
                        case renderMode_t::realGradient:value = dx.real(); break;
                        case renderMode_t::imagGradient:value = dx.imag(); break;
                        case renderMode_t::absGradient:value = std::abs(dx); break;
                        case renderMode_t::f2:value = f; break;
                        case renderMode_t::Jx:value = J.dfdx; break;
                        case renderMode_t::Jy:value = J.dfdy; break;
                        // case renderMode_t::realGradient:value = std::abs(fx); break;
                    }

                    scalarFieldData[y * nx + x] = (float) value;
        }
        else{

int32_t i = 0;
        cheb::complex pos = cheb::complex(xPos,yPos);
        cheb::complex prior = pos;
        inputFieldData[y * nx + x] = pos;
        // newton
        auto [state, location, positions, values, steps] = optimize(pos, method);
        pos = location;
        //auto val = std::abs(values[values.size() - 1]);
        auto fn = evalFunction(pos);
        auto iters = (int32_t) positions.size();

        auto cx = std::clamp(pos.real(), -1., 1.);
        auto cy = std::clamp(pos.imag(), -1., 1.);

        auto val = (float)std::sqrt(cx * cx + cy * cy);

        scalarFieldData[y * nx + x] = val;
        vectorFieldData[y * nx + x] = pos;
        functionFieldData[y*nx + x] = fn;
        iterationFieldData[y * nx + x] = state;

        // scalarFieldData[y * nx + x] = (float)std::abs(evalFunction(pos));
        if(clustering)
            scalarFieldData[y * nx + x] = std::abs(-evalFunction(pos) / evalDerivative(pos));
}
        }
        int32_t iter = ++iterationCounter;
        if(iter % 128 == 0){
            #pragma omp critical
            {
                pbl.update(iter,1);
            }
        }
    }
    std::cout << "\nFinished generating Visualization" << std::endl;

    auto min = *std::min_element(scalarFieldData, scalarFieldData + nx * ny);
    auto max = *std::max_element(scalarFieldData, scalarFieldData + nx * ny);

    auto minNorm = DBL_MAX;
    auto maxNorm = 0.;
    for (int32_t y = 0; y < screenHeight; ++y) {
        for (int32_t x = 0; x < screenWidth; ++x) {
            int32_t xi = std::clamp((int32_t)::floor(((scalar) x) / ((scalar) screenWidth) * (scalar) nx),0,nx - 1);
            int32_t yi = std::clamp((int32_t)::floor(((scalar) y) / ((scalar) screenHeight) * (scalar) ny),0,ny - 1);
            
            auto cur = vectorFieldData[yi * nx + xi];
            if(cur == cur && !std::isinf(cur.real())){
                minNorm = std::min(minNorm, std::abs(cur));
                maxNorm = std::max(maxNorm, std::abs(cur));                
            }
        }
    }
    maxNorm = 1.;
    for (int32_t y = 0; y < screenHeight; ++y) {
        for (int32_t x = 0; x < screenWidth; ++x) {
            int32_t xi = std::clamp((int32_t)::floor(((scalar) x) / ((scalar) screenWidth) * (scalar) nx),0,nx - 1);
            int32_t yi = std::clamp((int32_t)::floor(((scalar) y) / ((scalar) screenHeight) * (scalar) ny),0,ny - 1);
            
            auto cur = vectorFieldData[yi * nx + xi];
            if(cur == cur && !std::isinf(cur.real())){
                        //   radial[y * screenWidth + x] = (std::clamp(std::abs(cur), minNorm, maxNorm) - minNorm )/ (maxNorm - minNorm);
                          radial[y * screenWidth + x] = (std::clamp(std::abs(cur), minNorm, maxNorm))/ (maxNorm);
                          angular[y * screenWidth + x] = std::arg(cur) + std::numbers::pi;    
            }
        }
    }



#pragma omp parallel for
    for (int32_t y = 0; y < screenHeight; ++y) {
        for (int32_t x = 0; x < screenWidth; ++x) {
            int32_t xi = std::clamp((int32_t)::floor(((scalar) x) / ((scalar) screenWidth) * (scalar) nx),0,nx - 1);
            int32_t yi = std::clamp((int32_t)::floor(((scalar) y) / ((scalar) screenHeight) * (scalar) ny),0,ny - 1);
            
            auto v = scalarFieldData[yi * nx + xi];
            v = (v - min) / (max - min);
            v = std::clamp(v, 0.f, 1.f);
            data[y * screenWidth + x] = v;
            if(mode == renderMode_t::fractal && scalarFieldData[yi * nx + xi] < -0.5)
                data[y * screenWidth + x] = -1.;
            if(scalarFieldData[yi * nx + xi] != scalarFieldData[yi * nx + xi] ||
            vectorFieldData[yi * nx + xi] != vectorFieldData[yi * nx + xi]  )
                data[y * screenWidth + x] = -1.;

            if(iterationFieldData[yi * nx + xi] == optimizationState::diverged){
                data[y * screenWidth + x] = -1.;
            }
            if(cycles)
            if(iterationFieldData[yi * nx + xi] == optimizationState::unstable){
                data[y * screenWidth + x] = -5.;
            }
        }
    }
    return std::make_pair(min, max);

}

void renderField() {
    static GLuint textureID, radialTextureID, angularTextureID, vao, texID, angTexID, radTexID, program, minUni, maxUni, texUnit, fracuni;
    static bool once = true;
    static float* data = new float[screenWidth * screenHeight];
    static float* angularData = new float[screenWidth * screenHeight];
    static float* radialData = new float[screenWidth * screenHeight];
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
uniform sampler2D radialTexture;
uniform sampler2D angularTexture;
uniform sampler1D colorRamp;
uniform float min = 0.0;
uniform float max = 1.0;
uniform int fractal = 1;
const float PI = 3.1415926535897932384626433832795;

vec3 hsl2rgb( vec3 c )
{
    vec3 rgb = clamp( abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),6.0)-3.0)-1.0, 0.0, 1.0 );

    return c.z + c.y * (rgb-0.5)*(1.0-abs(2.0*c.z-1.0));
}

void main(){
float v = texture( renderedTexture, UV ).r;
float rel = v;//(v - min) / (max - min);
if(rel < -0.5){
    if(rel < -2.5)
        color = vec3(0,1,0);
    else
        color = vec3(1,0,0);
    }else{
        if(fractal == 0){
        float clamped = clamp(rel, 0.0, 1.0);
        vec3 col = texture(colorRamp, clamped).xyz;
        color = col;
    }else{
        float hue = texture(angularTexture, UV).r / (2.f * PI);
        float saturation = 1.;
        float lightness = texture(radialTexture,UV).r;
        lightness = lightness * 0.5f + 0.5f;
        vec3 hsl = vec3(hue, saturation, lightness);
        //vec3 col = hsl2rgb(hsl);
        vec3 col = texture(colorRamp, hue).xyz * lightness;
        //vec3 col = vec3(hue, hue, hue);
        color = col;
        if(texture(radialTexture,UV).r < 1e-3f)
        color = vec3(0,0,0);
        }

    }
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
    glGenTextures(1, &radialTextureID);
    glGenTextures(1, &angularTextureID);

    // "Bind" the newly created texture : all future texture functions will modify this texture
    glBindTexture(GL_TEXTURE_2D, textureID);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, screenWidth, screenHeight, 0, GL_RED, GL_FLOAT, data);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, radialTextureID);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, screenWidth, screenHeight, 0, GL_RED, GL_FLOAT, data);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, angularTextureID);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, screenWidth, screenHeight, 0, GL_RED, GL_FLOAT, data);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    GLuint quad_vertexbuffer;
    glGenBuffers(1, &quad_vertexbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, quad_vertexbuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_quad_vertex_buffer_data), g_quad_vertex_buffer_data, GL_STATIC_DRAW);

    // Create and compile our GLSL program from the shaders
    program = createProgram(vtxShader, frgShader);

    glUseProgram(program);
    angTexID = glGetUniformLocation(program, "angularTexture");
    radTexID = glGetUniformLocation(program, "radialTexture");
    texID = glGetUniformLocation(program, "renderedTexture");
    minUni = glGetUniformLocation(program, "min");
    maxUni = glGetUniformLocation(program, "max");
    fracuni = glGetUniformLocation(program, "fractal");
    //std::cout << "texID: " << texID << std::endl;
    //std::cout << "minUni: " << minUni << std::endl;
    //std::cout << "maxUni: " << maxUni << std::endl;

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, textureID);
    glUniform1i(texID, 1);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, radialTextureID);
    glUniform1i(radTexID, 2);
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, angularTextureID);
    glUniform1i(angTexID, 3);
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
    std::string file_name = "twilight.png";
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
    try{
    std::string file_name = resolveFile(std::string("cfg/") + ParameterManager::instance().get<std::string>("colorMap.colorMap") + ".png").string();

    if (std::filesystem::exists(file_name)) {
        std::cout << "Loading " << file_name << std::endl;
        unsigned char* image_data = stbi_load(file_name.c_str(), &image_width, &image_height, NULL, 4);
        delete[] img;
        img = new v4[image_width];
        for (int32_t it = 0; it < image_width; ++it) {
            img[it] = v4{
            (float)image_data[it * 4 + 0 + image_width * (image_height / 2)],
            (float)image_data[it * 4 + 1 + image_width * (image_height / 2)],
            (float)image_data[it * 4 + 2 + image_width * (image_height / 2)], 255.f };
            //std::cout << it << ": [ " << img[it].x << " " << img[it].y << " " << img[it].z <<
            //    " ]" << std::endl;
        }
        //img = QImage(QString::fromStdString(file_name));
        //img.load(QString(file_name.c_str()));
        //std::cout << image_width << " : " << image_height << std::endl;
    }}
    catch (...) {}
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
  WATCH(scalar, adam, alpha);
  WATCH(scalar, adam, beta1);
  WATCH(scalar, adam, beta2);
  WATCH(scalar, adam, eps);

  WATCH(scalar, domain, xmin);
  WATCH(scalar, domain, xmax);
  WATCH(scalar, domain, ymin);
  WATCH(scalar, domain, ymax);
    WATCH(std::string, field, method);
    WATCH(scalar, field, learningRate);
    WATCH(bool, field, clustering);
    WATCH(bool, field, cycles);
static int32_t oldIndex = -1;
auto index = ParameterManager::instance().get<int32_t>("index");
dirty = dirty || (index != oldIndex);
    if (dirty) {
        oldIndex = index;
        auto coefficients = globalCoefficients[index];
        //globalFunction = cheb::Function(std::vector<cheb::IntervalFunction>{cheb::IntervalFunction(coefficients)});

        auto [min,max] = updateField(data, angularData, radialData);
        //std::cout << buff << " -> " << vMode << std::endl;
        //std::cout << "New min and max: " << min << " x " << max << std::endl;
        ParameterManager::instance().get<scalar>("field.min") = min;
        ParameterManager::instance().get<scalar>("field.max") = max;
        glUseProgram(program);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, screenWidth, screenHeight, GL_RED, GL_FLOAT, data);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, radialTextureID);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, screenWidth, screenHeight, GL_RED, GL_FLOAT, radialData);
        glActiveTexture(GL_TEXTURE3);
        glBindTexture(GL_TEXTURE_2D, angularTextureID);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, screenWidth, screenHeight, GL_RED, GL_FLOAT, angularData);
        glUseProgram(0);
    }
    glUseProgram(program);
    //std::cout << textureID << " -> " << texUnit << std::endl;
    glBindVertexArray(vao);
    glUniform1i(fracuni, ParameterManager::instance().get<std::string>("colorMap.renderMode") == "fractal" ? 1 : 0);
    glUniform1f(minUni, ParameterManager::instance().get<scalar>("field.min"));
    glUniform1f(maxUni, ParameterManager::instance().get<scalar>("field.max"));
    //GLuint samplerLocation = glGetUniformLocation(program, "colorRamp");
    //glUniform1i(samplerLocation, 0);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, textureID);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, radialTextureID);
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, angularTextureID);
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
        static  bool once = true;
    bool dirty = once;
//     WATCH(double, path, x);
//     WATCH(double, path, y);
//     WATCH(std::string, field, method);
//   WATCH(scalar, adam, alpha);
//   WATCH(scalar, adam, beta1);
//   WATCH(scalar, adam, beta2);
//   WATCH(scalar, adam, eps);
//     WATCH(scalar, field, learningRate);
//     if(dirty){
//         trace.clear();

//         once = false;
// 		auto x =	ParameterManager::instance().get<scalar>("path.x");
// 		auto y =	ParameterManager::instance().get<scalar>("path.y");

// auto learningRate = std::pow(10., ParameterManager::instance().get<scalar>("field.learningRate"));
// auto stringMethod = ParameterManager::instance().get<std::string>("field.method");
// optimizationMethod method = getMethod(stringMethod);



//         cheb::complex pos = cheb::complex(x,y);
//         std::cout << "Starting path tracing at " << pos.real() << " + " << pos.imag() << "i\n";
//         cheb::complex prior = pos;
//         auto [state, location, positions, values, steps] = optimize(pos, method);
//         // auto [state, location, positions, values, steps] = BFGS(pos);
//         for(int32_t i = 0; i < positions.size(); ++ i){
//             trace.push_back(std::make_tuple(values[i], steps[i], positions[i]));
//             printf("\t[%03d]: f(%g + %gi) = %g + %gi -> %g + %gi\n", i, positions[i].real(), positions[i].imag(), values[i].real(), values[i].imag(), steps[i].real(), steps[i].imag());
//         }

//         std::cout << "Final path position : " << location.real() << " + " << location.imag() << std::endl;
//     }
    auto [domainMin, domainMax] = getDomain();
    glLoadIdentity();
    glUseProgram(0);

    glOrtho(domainMin.x(), domainMax.x(), domainMin.y(), domainMax.y(), 0, 1);
    // glOrtho(domainMin.x(), domainMax.x(), domainMax.y(), domainMin.y(), 0, 1);
    // glBegin(GL_LINES);
    // for(int32_t i = 0; i < trace.size() - 1; ++i){
    //     auto [fxl,dxl,pl] = trace[i];
    //     auto [fxr,dxr,pr] = trace[i+1];
    //     glVertex2f(pl.real(), pl.imag());
    //     glVertex2f(pr.real(), pr.imag());
    // }

    glEnd();


    glColor4f(0.8f, 0.f, 0.f, 1);
}