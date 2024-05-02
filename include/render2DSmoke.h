#ifndef RENDER_2D_SMOKE_H
#define RENDER_2D_SMOKE_H

#include <GLFW/glfw3.h>

namespace Fluid {

void render_2D_Smoke(GLFWwindow* window, double *density, int grid_width, int grid_height, int render_cell_size);

} // Fluid namespace

#endif // RENDER_2D_SMOKE_H