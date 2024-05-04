#include "viewer_2d.h"

#include <GLFW/glfw3.h>

Viewer2D::Viewer2D(int grid_width, int grid_height, int cell_size)
    : _grid_w(grid_width), _grid_h(grid_height), _cell_sz(cell_size)
{
}

void Viewer2D::render(std::vector<std::vector<double>> density) {
    // Orthographic projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, _grid_w*_cell_sz, _grid_h*_cell_sz, 0, -1, 1);
    // Render each cell as a quad
    glBegin(GL_QUADS);
    for (int x = 0; x < _grid_w; x ++) {
        for (int y = 0; y < _grid_h; y ++) {
            double d = density[x][y];
            if (d > 0.0) glColor3f(0.0, 0.0, d);
            else glColor3f(1.0, 1.0, 1.0);
            glVertex2f(x*_cell_sz, y*_cell_sz);
            glVertex2f((x+1)*_cell_sz, y*_cell_sz);
            glVertex2f((x+1)*_cell_sz, (y+1)*_cell_sz);
            glVertex2f(x*_cell_sz, (y+1)*_cell_sz);
        }
    }
    glEnd();
}