#include "render2DSmoke.h"
#include <GLFW/glfw3.h>

void render_2D_Smoke(GLFWwindow* window, double *density, int grid_width, int grid_height, int render_cell_size) {
    // Orthographic projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, grid_width * render_cell_size, grid_height * render_cell_size, 0, -1, 1);
    // Render each cell as a quad
    glBegin(GL_QUADS);
    for (int y = 0; y < grid_height; y ++) {
        for (int x = 0; x < grid_width; x ++) {
            double rho = density[y * grid_width + x];
            if (rho > 0.0) glColor3f(0.0, 0.0, rho);
            else glColor3f(1.0, 1.0, 1.0);
            glVertex2f(x*render_cell_size, y*render_cell_size);
            glVertex2f((x+1)*render_cell_size, y*render_cell_size);
            glVertex2f((x+1)*render_cell_size, (y+1)*render_cell_size);
            glVertex2f(x*render_cell_size, (y+1)*render_cell_size);
        }
    }
    glEnd();
}