#include "smoke_solver_2d.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "config_reader.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <stdio.h>
#define GL_SILENCE_DEPRECATION
#if defined(IMGUI_IMPL_OPENGL_ES2)
#include <GLES2/gl2.h>
#endif
#include <GLFW/glfw3.h> // Will drag system OpenGL headers

#if defined(_MSC_VER) && (_MSC_VER >= 1900) && !defined(IMGUI_DISABLE_WIN32_FUNCTIONS)
#pragma comment(lib, "legacy_stdio_definitions")
#endif

#include <iostream>

#include <time.h>
#include <stdlib.h>

static void glfw_error_callback(int error, const char* description)
{
    fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}

int main(int argc, char* argv[])
{
    // default config file
    std::string config_filename = "../config/example.yaml";
    if (argc > 1) { // user-sepcified config file 
        config_filename = argv[1];
    }
    SolverConfig solver_config = read_YAML(config_filename);

    int width = solver_config.width;
    int height = solver_config.height;

    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
        return 1;

    // Decide GL+GLSL versions
#if defined(IMGUI_IMPL_OPENGL_ES2)
    // GL ES 2.0 + GLSL 100
    const char* glsl_version = "#version 100";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_ES_API);
#elif defined(__APPLE__)
    // GL 3.2 + GLSL 150
    const char* glsl_version = "#version 150";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // Required on Mac
#else
    // GL 3.3 + GLSL 330
    const char* glsl_version = "#version 330";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // 3.0+ only
#endif


    // Create window with graphics context
    GLFWwindow* window = glfwCreateWindow(width, height, "Dear ImGui GLFW+OpenGL3 example", nullptr, nullptr);
    if (window == nullptr)
        return 1;
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsLight();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
#ifdef __EMSCRIPTEN__
    ImGui_ImplGlfw_InstallEmscriptenCanvasResizeCallback("#canvas");
#endif
    ImGui_ImplOpenGL3_Init(glsl_version);

    // Generating directory storing frames
    //  with current time as name
    time_t rawtime;
    struct tm * timeinfo;
    char dirname[256];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    sprintf(dirname, "./output/%d-%02d-%02d_%02d-%02d-%02d",
            timeinfo->tm_year + 1900, timeinfo->tm_mon + 1, timeinfo->tm_mday,
            timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
    // Creating the directory
    char mkdir_cmd[512];
    sprintf(mkdir_cmd, "mkdir -p %s", dirname);
    system(mkdir_cmd);

    // Mix the smoke color to the background (white) color
    auto mix_to_white = [](int component, double ratio) {
        int x = 100-std::min(static_cast<int>(ratio*100+0.5), 100);
        int rounded = std::min((((100-x)*component+x*255)/100), 255);
        int clipped = std::max(0, std::min(rounded, 255));
        return static_cast<unsigned char>(clipped);
    };

    // Smoke solver
    SmokeSolver2D smoke_solver(
        solver_config.width, solver_config.height, solver_config.dx,
        solver_config.alpha, solver_config.beta,
        solver_config.amb_T, solver_config.amb_s,
        solver_config.wind_u, solver_config.wind_v,
        solver_config.rate_T, solver_config.rate_s, solver_config.T_target,
        solver_config.int_sch
    );

    int frame = 0;

    // Main loop
#ifdef __EMSCRIPTEN__
    // For an Emscripten build we are disabling file-system access, so let's not attempt to do a fopen() of the imgui.ini file.
    // You may manually call LoadIniSettingsFromMemory() to load settings from your own storage.
    io.IniFilename = nullptr;
    EMSCRIPTEN_MAINLOOP_BEGIN
#else
    while (!glfwWindowShouldClose(window))
#endif
    {
        glfwPollEvents();

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // Parameter Window 
        {
            ImGui::Begin("Parameters");
            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
            ImGui::End();
        }

        // Simulate one step
        // TODO: dt
        smoke_solver.step();

        // Rendering
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        // Write this frame to .png
        unsigned char* image_data = (unsigned char*)malloc(width * height * 3);
        for (int y = 0; y < height; y ++) {
            for (int x = 0; x < width; x ++) {
                int index = ((height-1-y) * width + x) * 3;
                double d = smoke_solver._s[x][y];
                image_data[index] = mix_to_white(0, d);     // Red
                image_data[index + 1] = mix_to_white(102, d);   // Green
                image_data[index + 2] = mix_to_white(204, d);   // Blue
            }
        }
        char filename[256];
        sprintf(filename, "%d.png", frame);
        char filepath[512];
        sprintf(filepath, "%s/%s", dirname, filename);
        stbi_write_png(filepath, width, height, 3, image_data, width * 3);
        // Free memory
        free(image_data);
        // Next frame
        frame ++;

        glfwSwapBuffers(window);
    }
#ifdef __EMSCRIPTEN__
    EMSCRIPTEN_MAINLOOP_END;
#endif

    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
