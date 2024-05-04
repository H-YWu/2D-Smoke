#ifndef VIEWER_2D_H 
#define VIEWER_2D_H 

#include <vector>

class Viewer2D {
public:
    Viewer2D(int grid_width, int grid_height, int cell_size=1);

    void render(std::vector<std::vector<double>> density);
private:
    int _grid_w, _grid_h;
    int _cell_sz;
};

#endif // VIEWER_2D_H