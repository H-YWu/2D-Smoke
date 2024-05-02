#include "simulate2DSmoke.h"

namespace Fluid {

void simulate_2D_Smoke(double *density, int grid_width, int grid_height) {
    for (int y = 0; y < grid_height; ++y) {
        for (int x = 0; x < grid_width; ++x) {
            if (y == grid_height - 1|| x == grid_width - 1) density[y * grid_width + x] = 0.5;
            else density[y * grid_width + x] = 0.0;
        }
    }
}

} // Fluid namespace