cmake_minimum_required(VERSION 3.10)
project(2D-Smoke)

set(CMAKE_CXX_STANDARD 11)

# Submodules
set(GLFW_DIR external/glfw)
add_subdirectory(${GLFW_DIR})

# ImGui {*.cpp,*.h} directories
set(IMGUI_DIR external/imgui)
set(IMGUI_BACKENDS_DIR ${IMGUI_DIR}/backends)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -I${IMGUI_DIR} -I${IMGUI_BACKENDS_DIR} -g -Wall -Wformat")

# Libraries
set(LIBS "")
# OpenGL libraries for Linux
set(LINUX_GL_LIBS -lGL)
# Set OpenGL ES libraries (if needed)
#set(OPENGL_ES_LIBS -lGLESv2)
# Find OpenGL library
find_package(OpenGL REQUIRED)

# Platform specific 
if(UNIX AND NOT APPLE)
    message("Linux")
    list(APPEND LIBS ${LINUX_GL_LIBS})
    find_package(PkgConfig REQUIRED)
    pkg_search_module(GLFW REQUIRED glfw3)
    list(APPEND LIBS ${GLFW_STATIC_LIBRARIES})
    include_directories(${GLFW_INCLUDE_DIRS})
endif()
if(APPLE)
    message("Mac OS X")
    list(APPEND LIBS "-framework OpenGL" "-framework Cocoa" "-framework IOKit" "-framework CoreVideo")
    include_directories(/usr/local/include /opt/local/include /opt/homebrew/include)
    link_directories(/usr/local/lib /opt/local/lib /opt/homebrew/lib)
    find_package(glfw3 REQUIRED)
    list(APPEND LIBS glfw)
endif()
if(WIN32)
    message("MinGW")
    list(APPEND LIBS glfw3 gdi32 opengl32 imm32)
endif()

# Source files
file(GLOB SOURCES
    src/*.cpp
)
# ImGui
file(GLOB IMGUI_SOURCES
    ${IMGUI_DIR}/*.cpp
)
# ImGui backends
file(GLOB IMGUI_BACKEND_SOURCES
    ${IMGUI_BACKENDS_DIR}/imgui_impl_glfw.cpp
    ${IMGUI_BACKENDS_DIR}/imgui_impl_opengl3.cpp
)
list(APPEND SOURCES ${IMGUI_SOURCES} ${IMGUI_BACKEND_SOURCES})

# Executable target
set(EXECUTABLE_NAME 2D-Smoke)
add_executable(${EXECUTABLE_NAME} ${SOURCES})

# Include directories
target_include_directories(${EXECUTABLE_NAME} PRIVATE ${IMGUI_DIR})
target_include_directories(${EXECUTABLE_NAME} PRIVATE ${IMGUI_BACKENDS_DIR})
target_include_directories(${EXECUTABLE_NAME} PRIVATE ${GLFW_DIR}/include)

# Link libraries
target_link_libraries(${EXECUTABLE_NAME} PRIVATE ${LIBS})
