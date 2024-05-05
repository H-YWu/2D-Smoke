#ifndef VIEWER_2D_H 
#define VIEWER_2D_H 

#include <vector>

class Viewer2D {
public:
    Viewer2D(int grid_width, int grid_height);

    void render(std::vector<std::vector<double>> density);
private:
    int _grid_w, _grid_h;
};

#endif // VIEWER_2D_H