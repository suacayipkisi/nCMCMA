#define _USE_MATH_DEFINES
#include <cmath>

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>

#include "gui.h"
#include "../engine/simEngine.h"

#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#if __has_include(<Eigen/Dense>)
#include <Eigen/Dense>
#include <Eigen/Geometry>
#elif __has_include(<eigen3/Eigen/Dense>)
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#endif

const char* vertexShaderSource = R"(
#version 330 core
layout (location = 0) in vec3 aPos;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
void main() {
    gl_Position = projection * view * model * vec4(aPos, 1.0);
}
)";

const char* fragmentShaderSource = R"(
#version 330 core
out vec4 FragColor;
uniform vec3 objectColor;
void main() {
    FragColor = vec4(objectColor, 1.0); 
}
)";

Eigen::Matrix4f getPerspective(float fovY, float aspect, float zNear, float zFar) {
    Eigen::Matrix4f P = Eigen::Matrix4f::Zero();
    float tanHalfFovy = std::tan(fovY / 2.0f);
    P(0, 0) = 1.0f / (aspect * tanHalfFovy);
    P(1, 1) = 1.0f / tanHalfFovy;
    P(2, 2) = -(zFar + zNear) / (zFar - zNear);
    P(3, 2) = -1.0f;
    P(2, 3) = -(2.0f * zFar * zNear) / (zFar - zNear);
    return P;
}

Eigen::Matrix4f getLookAt(const Eigen::Vector3f& eye, const Eigen::Vector3f& center, const Eigen::Vector3f& up) {
    Eigen::Vector3f f = (center - eye).normalized();
    Eigen::Vector3f u = up.normalized();
    Eigen::Vector3f s = f.cross(u).normalized();
    u = s.cross(f);
    Eigen::Matrix4f res = Eigen::Matrix4f::Identity();
    res(0, 0) = s.x(); res(0, 1) = s.y(); res(0, 2) = s.z();
    res(1, 0) = u.x(); res(1, 1) = u.y(); res(1, 2) = u.z();
    res(2, 0) = -f.x(); res(2, 1) = -f.y(); res(2, 2) = -f.z();
    res(0, 3) = -s.dot(eye);
    res(1, 3) = -u.dot(eye);
    res(2, 3) = f.dot(eye);
    return res;
}

unsigned int compileShader(unsigned int type, const char* source) {
    unsigned int id = glCreateShader(type);
    glShaderSource(id, 1, &source, nullptr);
    glCompileShader(id);
    int result;
    glGetShaderiv(id, GL_COMPILE_STATUS, &result);
    if (!result) {
        int logLength = 0;
        glGetShaderiv(id, GL_INFO_LOG_LENGTH, &logLength);
        std::vector<char> infoLog(logLength);
        glGetShaderInfoLog(id, logLength, &logLength, infoLog.data());
        std::cerr << "Shader Hatasi: " << infoLog.data() << '\n';
        glDeleteShader(id);
        return 0;
    }
    return id;
}

void initAndRunGui() {
    if (!glfwInit()) {
        std::cerr << "GLFW coudn't started!\n";
        return;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(1280, 720, "nCMCMA Modal Analysis", nullptr, nullptr);
    if (!window) { 
        std::cerr << "Window couldn't created!\n";
        glfwTerminate(); 
        return; 
    }
    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "GLAD couldn't loaded!\n";
        glfwDestroyWindow(window);
        glfwTerminate();
        return;
    }

    glfwSwapInterval(1); 

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330 core");

    ImGuiIO& io = ImGui::GetIO();
    io.FontGlobalScale = 1.15f; 

    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = 0.0f;
    style.FrameRounding = 6.0f; 
    style.GrabRounding = 6.0f;  
    style.FramePadding = ImVec2(8.0f, 6.0f); 
    style.ItemSpacing = ImVec2(8.0f, 10.0f); 
    style.Colors[ImGuiCol_WindowBg] = ImVec4(0.12f, 0.12f, 0.14f, 1.0f);
    style.Colors[ImGuiCol_Header] = ImVec4(0.2f, 0.2f, 0.25f, 1.0f);

    unsigned int program = glCreateProgram();
    unsigned int vs = compileShader(GL_VERTEX_SHADER, vertexShaderSource);
    unsigned int fs = compileShader(GL_FRAGMENT_SHADER, fragmentShaderSource);
    glAttachShader(program, vs);
    glAttachShader(program, fs);
    glLinkProgram(program);

    int linkStatus;
    glGetProgramiv(program, GL_LINK_STATUS, &linkStatus);
    if (!linkStatus) {
        int logLen = 0;
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLen);
        std::vector<char> log(logLen);
        glGetProgramInfoLog(program, logLen, nullptr, log.data());
        std::cerr << "Program Link Error: " << log.data() << '\n';
        return;
    }

    glDeleteShader(vs);
    glDeleteShader(fs);

    std::vector<float> vertices;
    std::vector<unsigned int> indices;
    const int sectors = 16;
    const int stacks = 16;
    for(int i = 0; i <= stacks; ++i) {
        float V = i / (float)stacks;
        float phi = V * (float)M_PI;
        for(int j = 0; j <= sectors; ++j) {
            float U = j / (float)sectors;
            float theta = U * ((float)M_PI * 2.0f);
            float x = std::cos(theta) * std::sin(phi);
            float y = std::cos(phi);
            float z = std::sin(theta) * std::sin(phi);
            vertices.push_back(x); vertices.push_back(y); vertices.push_back(z);
        }
    }
    for(int i = 0; i < stacks; ++i) {
        for(int j = 0; j < sectors; ++j) {
            indices.push_back(i * (sectors + 1) + j);
            indices.push_back((i + 1) * (sectors + 1) + j);
            indices.push_back((i + 1) * (sectors + 1) + j + 1);
            indices.push_back(i * (sectors + 1) + j);
            indices.push_back((i + 1) * (sectors + 1) + j + 1);
            indices.push_back(i * (sectors + 1) + j + 1);
        }
    }

    unsigned int VAO, VBO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glEnable(GL_DEPTH_TEST);

    int modelLoc = glGetUniformLocation(program, "model");
    int viewLoc = glGetUniformLocation(program, "view");
    int projLoc = glGetUniformLocation(program, "projection");
    int colorLoc = glGetUniformLocation(program, "objectColor");

    float animScale = 1.0f; 
    static int selectedModeIndex = 0;

    if (g_params.massNum <= 0) g_params.massNum = 3;
    if (g_params.radius <= 0.0f) g_params.radius = 0.05f;

    float camRadius = g_params.massNum * g_params.radius * 6.0f;
    float camYaw = (float)M_PI / 4.0f;
    float camPitch = (float)M_PI / 6.0f;
    Eigen::Vector3f camTarget(0.0f, 0.0f, 0.0f);

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        const float panelWidth = 500.0f; 

        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(panelWidth, (float)display_h));
        
        ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | 
                                        ImGuiWindowFlags_NoCollapse | 
                                        ImGuiWindowFlags_NoResize | 
                                        ImGuiWindowFlags_NoMove | 
                                        ImGuiWindowFlags_NoBringToFrontOnFocus;

        ImGui::Begin("Control Panel", nullptr, window_flags);
        
        ImGui::TextColored(ImVec4(0.4f, 0.8f, 1.0f, 1.0f), "nCMCMA MODAL ANALYSIS");
        ImGui::Separator();
        ImGui::Spacing();

        if (ImGui::CollapsingHeader("Lattice and Material Constants", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::InputInt("Lattice Size (N)", &g_params.massNum);
            ImGui::InputFloat("Sphere Radius (m)", &g_params.radius, 0.001f, 0.01f, "%.4f");
            ImGui::InputFloat("Density (kg/m^3)", &g_params.density);
            ImGui::InputFloat("Spring Constant (N/m)", &g_params.stiffnessConst);
        }

        if (ImGui::CollapsingHeader("Rayleigh Damping Ratios", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::InputFloat("Alpha (Mass)", &g_params.alpha, 0.001f, 0.01f, "%.5f");
            ImGui::InputFloat("Beta (Spring)", &g_params.beta, 0.000001f, 0.00001f, "%.8f");
        }

        if (ImGui::CollapsingHeader("Simulation Control", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::SliderFloat("Visual Motion Multiplier", &animScale, 0.1f, 100.0f);
            ImGui::SliderFloat("Max Frequency (Hz)", &g_params.maxTargetHz, 1.0f, 500.0f);

            if (g_simState.isRunning) {
                ImGui::ProgressBar(g_simState.progress.load(), ImVec2(0.0f, 0.0f));
                ImGui::BeginDisabled(); ImGui::Button("Calculating..."); ImGui::EndDisabled();
            } else {
                if (ImGui::Button("Start Simulation")) {
                    selectedModeIndex = 0;
                    std::thread(runSimulation).detach();
                }
            }
        }

        std::vector<double> availableFrequencies;
        {
            std::lock_guard<std::mutex> lock(g_simState.dataMutex);
            for (const auto& res : g_simState.results) {
                availableFrequencies.push_back(res.hz);
            }
        }

        if (!availableFrequencies.empty()) {
            if (selectedModeIndex >= (int)availableFrequencies.size()) {
                selectedModeIndex = (int)availableFrequencies.size() - 1;
            }
            
            std::string previewValue = std::to_string(availableFrequencies[selectedModeIndex]) + " Hz";
            if (ImGui::BeginCombo("Visualized Mode", previewValue.c_str())) {
                for (size_t n = 0; n < availableFrequencies.size(); n++) {
                    const bool is_selected = (selectedModeIndex == (int)n);
                    std::string itemLabel = "Mode " + std::to_string(n+1) + " : " + std::to_string(availableFrequencies[n]) + " Hz";
                    if (ImGui::Selectable(itemLabel.c_str(), is_selected)) {
                        selectedModeIndex = (int)n;
                    }
                    if (is_selected) ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }
        } else {
            ImGui::Text("No calculated modes yet.");
        }

        ImGui::End();

        if (!io.WantCaptureMouse) { 
            if (ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
                camYaw -= io.MouseDelta.x * 0.005f;
                camPitch += io.MouseDelta.y * 0.005f;
                if (camPitch > (float)M_PI / 2.0f - 0.01f) camPitch = (float)M_PI / 2.0f - 0.01f;
                if (camPitch < -(float)M_PI / 2.0f + 0.01f) camPitch = -(float)M_PI / 2.0f + 0.01f;
            }
            camRadius -= io.MouseWheel * g_params.radius * 2.0f;
            if (camRadius < g_params.radius * 2.0f) camRadius = g_params.radius * 2.0f;
        }

        glViewport(0, 0, display_w, display_h);
        glClearColor(0.08f, 0.08f, 0.08f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        int renderWidth = display_w - (int)panelWidth;
        if (renderWidth < 1) renderWidth = 1;
        glViewport((int)panelWidth, 0, renderWidth, display_h);

        float camX = camRadius * std::cos(camPitch) * std::cos(camYaw);
        float camY = camRadius * std::sin(camPitch);
        float camZ = camRadius * std::cos(camPitch) * std::sin(camYaw);
        Eigen::Vector3f eye(camX, camY, camZ);

        float fovRad = 45.0f * ((float)M_PI / 180.0f);
        Eigen::Matrix4f projection = getPerspective(fovRad, (float)renderWidth / (float)display_h, 0.1f, 100.0f);
        Eigen::Matrix4f view = getLookAt(eye, camTarget, Eigen::Vector3f(0, 1, 0));
        
        glUseProgram(program);
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, view.data());
        glUniformMatrix4fv(projLoc, 1, GL_FALSE, projection.data());

        glBindVertexArray(VAO);

        float t = (float)glfwGetTime() * 5.0f;
        float spacing = g_params.radius * 3.0f;
        int N = g_params.massNum;
        int N2 = N * N;

        std::vector<std::complex<double>> currentDisp;
        {
            std::lock_guard<std::mutex> lock(g_simState.dataMutex);
            if (selectedModeIndex < (int)g_simState.results.size()) {
                currentDisp = g_simState.results[selectedModeIndex].displacements;
            }
        }

        bool hasData = (currentDisp.size() == (size_t)(6 * N * N * N));

        float maxDispMag = 0.00001f;
        if (hasData) {
            for (int i = 0; i < N * N * N; ++i) {
                int baseIdx = i * 6;
                double dx = currentDisp[baseIdx + 0].real();
                double dy = currentDisp[baseIdx + 1].real();
                double dz = currentDisp[baseIdx + 2].real();
                float mag = (float)std::sqrt(dx * dx + dy * dy + dz * dz);
                if (mag > maxDispMag) {
                    maxDispMag = mag;
                }
            }
        }
        
        for (int i = 0; i < N * N * N; ++i) {
            int ix = i % N;
            int iy = (i % N2) / N;
            int iz = i / N2;

            Eigen::Vector3f basePos(
                (ix - N/2.0f) * spacing, 
                (iy - N/2.0f) * spacing, 
                (iz - N/2.0f) * spacing
            );

            Eigen::Vector3f dynamicTrans(0, 0, 0);
            Eigen::Vector3f dynamicRot(0, 0, 0);
            float dispMag = 0.0f;

            if (hasData) {
                int baseIdx = i * 6;
                for(int d = 0; d < 3; ++d) {
                    auto disp = currentDisp[baseIdx + d];
                    dynamicTrans[d] = (float)(disp.real() * std::cos(t) - disp.imag() * std::sin(t)) * animScale;
                }
                for(int d = 0; d < 3; ++d) {
                    auto rot = currentDisp[baseIdx + 3 + d];
                    dynamicRot[d] = (float)(rot.real() * std::cos(t) - rot.imag() * std::sin(t)) * animScale;
                }
                double dx = currentDisp[baseIdx + 0].real();
                double dy = currentDisp[baseIdx + 1].real();
                double dz = currentDisp[baseIdx + 2].real();
                dispMag = (float)std::sqrt(dx * dx + dy * dy + dz * dz);
            }

            Eigen::Affine3f model = Eigen::Affine3f::Identity();
            model.translate(basePos + dynamicTrans);
            
            float rotNorm = dynamicRot.norm();
            if (rotNorm > 1e-6f) {
                model.rotate(Eigen::AngleAxisf(rotNorm, dynamicRot.normalized()));
            }
            model.scale(g_params.radius);

            float r, g, b;
            if (hasData) {
                float normT = std::clamp(dispMag / maxDispMag, 0.0f, 1.0f);
                r = std::clamp(2.0f * normT - 0.5f, 0.0f, 1.0f);
                g = std::clamp(1.0f - std::abs(2.0f * normT - 1.0f), 0.0f, 1.0f);
                b = std::clamp(1.5f - 2.0f * normT, 0.0f, 1.0f);
            } else {
                r = (N > 1) ? static_cast<float>(ix) / (N - 1) : 1.0f;
                g = (N > 1) ? static_cast<float>(iy) / (N - 1) : 1.0f;
                b = (N > 1) ? static_cast<float>(iz) / (N - 1) : 1.0f;
                r = 0.2f + 0.8f * r;
                g = 0.2f + 0.8f * g;
                b = 0.2f + 0.8f * b;
            }
            
            glUniform3f(colorLoc, r, g, b);
            glUniformMatrix4fv(modelLoc, 1, GL_FALSE, model.matrix().data());
            glDrawElements(GL_TRIANGLES, (GLsizei)indices.size(), GL_UNSIGNED_INT, 0);
        }

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
}