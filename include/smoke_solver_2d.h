#ifndef SMOKE_SOLVER_2D_H
#define SMOKE_SOLVER_2D_H

#include <vector>
#include <utility>

enum CellLabel {SOLID, FLUID, EMPTY};

enum FieldType {U, V, ST, SS};

void simulate_2D_Smoke(double *density, int width, int height);

// Simulate smoke in the smoke-free open air
//  with constant air velocity
class SmokeSolver2D {
public:
    // Simulation grid: width x height
    // Cell size: dx x dx
    // Smoke bouyance parameters: alpha, beta
    // Air parameters: ambient tempreature and ambient smoke concentration
    // Open air velocity: (wind_u, wind_v)
    // Smoke emitter parameters: rate of increase of tempreature and smoke concentration, target tempreature
    SmokeSolver2D(
        int grid_width, int grid_height, double dx,
        double alpha=0.5, double beta=0.1,
        double ambient_T=273.0, double ambient_s=0.0,
        double wind_u=1.0, double wind_v=0.0,
        double rate_T=8.0, double rate_s=10.0, double T_target=300.0
    );

    // Simulate smoke by one step
    void step();

private:
    // Simulation setup
    void init();

    // Calculate max dt from CFL condition
    void CFL_dt(double &dt);
    // Add smoke sources by one step
    void step_source(double dt);

    // Advect the field 
    void advect(FieldType ft, double dt);
    // Add external forces, esp. gravity and bouyancy
    void force(double dt);
    // Ensuring incompressible condition
    void project(double dt);

    void backward_Euler(FieldType ft, int x_G, int y_G, double dt, double &x_P, double &y_P);
    void RK2(FieldType ft, int x_G, int y_G, double dt, double &x_P, double &y_P);

    void solve_pressure(double dt);
    void update_uv_incompressible(double dt);

    // Interpolate velocity (u, v) at an arbitrary grid point (x, y)
    void blerp_uv(double x, double y, double &u, double &v);
    // Interpolate velocity (u, v) at a regular grid point (x, y)
    void cell_uv(FieldType ft, int x, int y, double &u, double &v);

    double u_with_bnd(int x, int y);
    double v_with_bnd(int x, int y);
    double T_with_bnd(int x, int y);
    double s_with_bnd(int x, int y);
    double p_with_bnd(int x, int y);

    CellLabel label(int x, int y);

private:
    static const double _g;
    static const double _rho0;  // TODO: variable rho
    // Grid resolution
    int _nx, _ny;
    // Cell size 
    double _dx;
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
public: // for rendering
    std::vector<std::vector<double>> _s;    // concentration
private:
    std::vector<std::vector<double>> _T_nxt;    // next step temperature
    std::vector<std::vector<double>> _s_nxt;    // next step concentration
    // For CFL 
    double _s_max, _dT_max, _u_max, _v_max; 
};

#endif // SIMULATE_2D_SMOKE_H