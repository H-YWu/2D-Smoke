#ifndef SIMULATE_2D_SMOKE_H
#define SIMULATE_2D_SMOKE_H

#include <vector>
#include <utility>

enum CellLabel {SOLID, FLUID, EMPTY};

enum QuantityType {U, V, ST, SS};

void simulate_2D_Smoke(double *density, int width, int height);

// simulating smoke in the open air we could assume
// a constant wind velocity U (perhaps zero) outside of the simulation domain.

class SmokeSolver2D {
public:
    SmokeSolver2D(
        int grid_width, int grid_height,
        int dx, int render_cell_size=1,
        double alpha=0.5, double beta=0.5,
        double ambient_T=273.0, double ambient_s=0.0,
        double wind_u=1.0, double wind_v=0.0,
        double rate_T=0.1, double rate_s=0.1, double T_target=300
    );

    void init();
    void step();
    void render();

private:
    void step_source(double dt);

    void advect(QuantityType qt, double dt, std::vector<std::vector<double>> &q_nxt);
    void force(double dt);
    void project(double dt);

    void backward_Euler(QuantityType qt, int x_G, int y_G, double dt, double &x_P, double &y_P);
    void RK2(QuantityType qt, int x_G, int y_G, double dt, double &x_P, double &y_P);

    void solve_pressure(double dt);
    void update_uv_incompressible(double dt);

    void blerp_uv(double x, double y, double &u, double &v);
    void cell_uv(QuantityType qt, int x_G, int y_G, double &u_G, double &v_G);
    double u_with_bnd(int x, int y);
    double v_with_bnd(int x, int y);
    double T_with_bnd(int x, int y);
    double s_with_bnd(int x, int y);
    double p_with_bnd(int x, int y);

    CellLabel label(int x, int y);

private:
    static const double _g;
    static const double _rho0;
    // Grid resolution
    int _nx, _ny;
    // Cell size 
    int _dx;
    // Rendering
    int _rc_sz;
    // Smoke parameters
    double _alpha, _beta;
    // Boundary
    double _amb_T, _amb_s;
    double _wind_u, _wind_v;
    // Smoke Sources
    double _r_T, _r_s, _T_tar;
    std::vector<std::pair<int,int>> _sources;
    // MAC cell
    std::vector<std::vector<CellLabel>> _label;
    std::vector<std::vector<double>> _pressure;
    std::vector<std::vector<double>> _u;    // velocity x
    std::vector<std::vector<double>> _v;    // velocity y
    std::vector<std::vector<double>> _u_nxt;    // next step velocity x
    std::vector<std::vector<double>> _v_nxt;    // next step velocity y
    std::vector<std::vector<double>> _T;    // temperature
    std::vector<std::vector<double>> _s;    // concentration
    std::vector<std::vector<double>> _T_nxt;    // next step temperature
    std::vector<std::vector<double>> _s_nxt;    // next step concentration
};

#endif // SIMULATE_2D_SMOKE_H