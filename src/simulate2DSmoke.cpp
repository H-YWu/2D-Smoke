#include "simulate2DSmoke.h"
#include "utils.h"
#include "pcg_solver.h"

#include <cmath>

void simulate_2D_Smoke(double *density, int grid_width, int grid_height) {
    for (int y = 0; y < grid_height; ++y) {
        for (int x = 0; x < grid_width; ++x) {
            if (y == grid_height - 1|| x == grid_width - 1) density[y * grid_width + x] = 0.5;
            else density[y * grid_width + x] = 0.0;
        }
    }
}

const double SmokeSolver2D::_g = -9.80665;   // m/s^2
const double SmokeSolver2D::_rho0 = 1.293;   // smoke-free air density, kg/m^{-3}

SmokeSolver2D::SmokeSolver2D(
    int grid_width, int grid_height,
    int dx, int render_cell_size,
    double alpha, double beta,
    // double ambient_T=273.0, double ambient_s=0.0,
    // double wind_u=1.0, double wind_v=0.0,
    // double rate_T=0.1, double rate_s=0.1, double T_target=300
    double ambient_T, double ambient_s,
    double wind_u, double wind_v,
    double rate_T, double rate_s, double T_target
)
    : _nx(grid_width), _ny(grid_height)
    , _dx(dx), _rc_sz(render_cell_size)
    , _alpha(alpha), _beta(beta)
    , _amb_T(ambient_T), _amb_s(ambient_s)
    , _wind_u(wind_u), _wind_v(wind_v)
    , _r_T(rate_T), _r_s(rate_s), _T_tar(T_target)
    , _label(grid_width, std::vector<CellLabel>(grid_height))
    , _pressure(grid_width, std::vector<double>(grid_height))
    , _u(grid_width+1, std::vector<double>(grid_height))
    , _v(grid_width, std::vector<double>(grid_height+1))
    , _u_nxt(grid_width+1, std::vector<double>(grid_height))
    , _v_nxt(grid_width, std::vector<double>(grid_height+1))
    , _T(grid_width, std::vector<double>(grid_height))
    , _s(grid_width, std::vector<double>(grid_height))
{
    init();
}

void SmokeSolver2D::init() {
    // Smoke sources
    _sources.emplace_back(std::make_pair(0, static_cast<int>(0.25 * _ny)));

    for (int x = 0; x <= _nx; x ++) {
        for (int y = 0; y <= _ny; y ++) {
            _label[x][y] = FLUID;
            _T[x][y] = 0.0;
            _s[x][y] = 0.0;
            _pressure[x][y] = 0.0;
        }
    }
    for (int x = 0; x <= _nx; x ++) {
        for (int y = 0; y < _ny; y ++) {
            _u[x][y] = 1.0;
            _u_nxt[x][y] = 1.0;
        }
    }
    for (int x = 0; x < _nx; x ++) {
        for (int y = 0; y <= _ny; y ++) {
            _v[x][y] = 0.0;
            _v_nxt[x][y] = 0.0;
        }
    }
}

void SmokeSolver2D::step() {
    // TODO: CFL dt
    double dt = 0.01;

    advect(U, dt, _u_nxt);
    advect(V, dt, _v_nxt);
    _u.swap(_u_nxt);
    _v.swap(_v_nxt);

    // diffuse();
    force(dt);
    project(dt);

    advect(ST, dt, _T_nxt);
    advect(SS, dt, _s_nxt);
    _T.swap(_T_nxt);
    _s.swap(_s_nxt);

    step_source(dt);
}


void SmokeSolver2D::step_source(double dt) {
    for (long unsigned int i = 0; i < _sources.size(); i ++) {
        int x = _sources[i].first;
        int y = _sources[i].second;
        _T[x][y] = _T[x][y] + (1.0 - std::exp(-_r_T * dt)) * (_T_tar - _T[x][y]);
        _s[x][y] = _s[x][y] + _r_s * dt;
    }
}

void SmokeSolver2D::advect(QuantityType qt, double dt, std::vector<std::vector<double>> &q_nxt) {
    int nx = _nx, ny = _ny;
    if (qt == U) nx ++;
    if (qt == V) ny ++;

    for (int x = 0; x < nx; x ++) {
        for (int y = 0; y < ny; y ++) {
            double x_P, y_P;
            // Default: Backward Euler
            // TODO: integration method select
            backward_Euler(qt, x, y, dt, x_P, y_P);
            // TODO: interpolation method select
            if (qt == U || qt == V) {
                double u_P, v_P;
                blerp_uv(x_P, y_P, u_P, v_P);
                if (qt == U) q_nxt[x][y] = u_P;
                if (qt == V) q_nxt[x][y] = v_P;
            }
            else {
                int x_low = static_cast<int>(std::floor(x_P));
                int y_low = static_cast<int>(std::floor(y_P));
                // Boundary condition
                if (x_low < 0 || x_low >= _nx || y_low < 0 || y_low >= _ny) {
                    if (qt == ST)
                        q_nxt[x][y] = T_with_bnd(x_low, y_low);
                    if (qt == SS)
                        q_nxt[x][y] = s_with_bnd(x_low, y_low);
                    continue;
                }
                double tx = (x_P-static_cast<double>(x_low))/1.0;
                double ty = (y_P-static_cast<double>(y_low))/1.0;
                if (qt == ST)
                    q_nxt[x][y] = blerp(
                        _T[x_low][y_low], _T[x_low+1][y_low],
                        _T[x_low+1][y_low], _T[x_low+1][y_low+1],
                        tx, ty);
                if (qt == SS)
                    q_nxt[x][y] = blerp(
                        _s[x_low][y_low], _s[x_low+1][y_low],
                        _s[x_low+1][y_low], _s[x_low+1][y_low+1],
                        tx, ty);
            }
        }
    }
}

void SmokeSolver2D::force(double dt) {
    // Gravity/Buoyancy only applies to v
    for (int x = 0; x < _nx; x ++) {
        for (int y = 0; y <= _ny; y ++) {
            double s = 0.5 * (s_with_bnd(x, y-1) + s_with_bnd(x, y));
            double T = 0.5 * (T_with_bnd(x, y-1) + T_with_bnd(x, y));
            double b = 
                (_alpha * s - _beta * (T - _amb_T)) * _g;
            _v_nxt[x][y] = _v[x][y] + dt * b;
        }
    }
    _v.swap(_v_nxt);
}

void SmokeSolver2D::project(double dt) {
    solve_pressure(dt); 
    // Update velocity of fluid from solved pressure
    update_uv_incompressible(dt);
}

void SmokeSolver2D::backward_Euler(QuantityType qt, int x_G, int y_G, double dt, double &x_P, double &y_P) {
    double u_G, v_G;
    cell_uv(qt, x_G, y_G, u_G, v_G);

    x_P = x_G - dt * u_G;
    y_P = y_G - dt * v_G;
}

void SmokeSolver2D::RK2(QuantityType qt, int x_G, int y_G, double dt, double &x_P, double &y_P) {
    double u_G, v_G;
    cell_uv(qt, x_G, y_G, u_G, v_G);

    double x_mid = x_G - 0.5 * dt * u_G;
    double y_mid = y_G - 0.5 * dt * v_G;

    double u_mid, v_mid;
    blerp_uv(x_mid, y_mid, u_mid, v_mid);

    x_P = x_G - dt * u_mid;
    y_P = y_G - dt * v_mid;
}

void SmokeSolver2D::solve_pressure(double dt) {
    // Construct the linear system
    int sz = _nx*_ny;
    SparseMatrix<double> A(sz);
    std::vector<double> result(sz);
    std::vector<double> b(sz);
    double scale_lhs = dt / (_rho0 * _dx * _dx);
    double scale_rhs = 1.0 / _dx;
    for (int x = 0; x < _nx; x ++) {
        for (int y = 0; y < _ny; y ++) {
            int r = x*_nx+y;
            if (label(x,y)==FLUID) {
                b[r] = -scale_rhs * (
                    _u[x+1][y] - _u[x][y] +
                    _v[x][y+1] - _v[x][y]
                );
                // i = (x, y)
                    // j = (x, y)
                A.add_to_element(r, r, 4.0 * scale_lhs);
                    // j = (x-1, y)
                if (label(x-1,y)==SOLID) {
                    A.add_to_element(r, r, -1.0 * scale_lhs);
                    // u(x,y) - usolid(x, y)
                    // Assume solids are all still
                    b[r] -= scale_rhs * (_u[x][y] - 0.0);
                } 
                else if (label(x-1,y)==FLUID) {
                    A.add_to_element(r, r-_nx, -1.0 * scale_lhs);
                } // EMPTY just 0
                    // j = (x+1, y)
                if (label(x+1,y)==SOLID) {
                    A.add_to_element(r, r, -1.0 * scale_lhs);
                    // u(x+1,y) - usolid(x+1, y)
                    // Assume solids are all still
                    b[r] += scale_rhs * (_u[x+1][y] - 0.0);
                } 
                else if (label(x+1,y)==FLUID) {
                    A.add_to_element(r, r+_nx, -1.0 * scale_lhs);
                } // EMPTY just 0
                    // j = (x, y-1)
                if (label(x,y-1)==SOLID) {
                    A.add_to_element(r, r, -1.0 * scale_lhs);
                    // v(x,y) - vsolid(x, y)
                    // Assume solids are all still
                    b[r] -= scale_rhs * (_v[x][y] - 0.0);
                }
                else if (label(x,y-1)==FLUID) {
                    A.add_to_element(r, r-1, -1.0 * scale_lhs);
                } // EMPTY just 0
                    // j = (x, y+1)
                if (label(x,y+1)==SOLID) {
                    A.add_to_element(r, r, -1.0 * scale_lhs);
                    // v(x,y+1) - vsolid(x, y+1)
                    // Assume solids are all still
                    b[r] += scale_rhs * (_v[x][y+1] - 0.0);
                } 
                else if (label(x,y+1)==FLUID) {
                    A.add_to_element(r, r+1, -1.0 * scale_lhs);
                } // EMPTY just 0
            } else {
                b[r] = 0.0;
            }
        }
    }
    // Solve the linear system
    PCGSolver<double> solver;
    double res;
    int iters;
    solver.solve(A, b, result, res, iters); 
    for (int x = 0; x < _nx; x ++) {
        for (int y = 0; y < _ny; y ++) {
            _pressure[x][y] = result[x*_nx+y];
        }
    }
    // TODO: check solver
}

void SmokeSolver2D::update_uv_incompressible(double dt) {
    double scale = dt / (_rho0 * _dx);
    for (int x = 0; x <= _nx; x ++) {
        for (int y = 0; y <= _ny; y ++) {
 
            if (y != _ny) { // update u
                if (label(x-1,y)==FLUID || label(x,y)==FLUID) {
                    if (label(x-1,y)==SOLID || label(x,y)==SOLID) {
                        // u(x,y) = usolid(x,y);
                        // Assume solids are all still
                        _u[x][y] = 0.0;
                    } else {
                        // TODO: pressure bnd
                        _u[x][y] -= scale * 
                            (_pressure[x][y] - _pressure[x-1][y]);
                    }
                } else {
                    // u(x,y) is unknown
                }
            }

            if (x != _nx) { // update v
                if (label(x,y-1)==FLUID || label(x,y)==FLUID) {
                    if (label(x,y-1)==SOLID || label(x,y)==SOLID) {
                        // v(x,y) = vsolid(x,y);
                        // Assume solids are all still
                        _v[x][y] = 0.0;
                    } else {
                        // TODO: pressure bnd
                        _v[x][y] -= scale * 
                            (_pressure[x][y] - _pressure[x][y-1]);
                    }
                } else {
                    // v(x,y) is unknown
                }
            }

        }
    }
}

void SmokeSolver2D::blerp_uv(double x, double y, double &u, double &v) {
    double xc = static_cast<double>(static_cast<int>(x+0.5));
    double yc = static_cast<double>(static_cast<int>(y+0.5));
    double u_cell[2][2], v_cell[2][2];
    for (int i = 0; i < 2; i ++) {
        for (int j = 0; j < 2; j ++) {
            u_cell[i][j] = 0.5*(u_with_bnd(xc+i,yc-1+j)+u_with_bnd(xc+i,yc+j));
            v_cell[i][j] = 0.5*(v_with_bnd(xc-1+i,yc+j)+v_with_bnd(xc+i,yc+j));
        }
    }
    double a = (x-xc+0.5)/1.0;
    double b = (y-yc+0.5)/1.0;
    u = blerp(u_cell[0][0], u_cell[1][0], u_cell[0][1], u_cell[1][1], a, b);
    v = blerp(v_cell[0][0], v_cell[1][0], v_cell[0][1], v_cell[1][1], a, b);
}

void SmokeSolver2D::cell_uv(QuantityType qt, int x_G, int y_G, double &u_G, double &v_G) {
    if (qt == U) {
        u_G = _u[x_G][y_G];
        v_G = 0.25 * (
            v_with_bnd(x_G-1, y_G) + v_with_bnd(x_G-1, y_G+1) +
            v_with_bnd(x_G, y_G) + v_with_bnd(x_G, y_G+1)
        );
    }
    if (qt == V) {
        u_G = 0.25 * (
            u_with_bnd(x_G-1, y_G) + u_with_bnd(x_G-1, y_G+1) +
            u_with_bnd(x_G, y_G) + u_with_bnd(x_G, y_G+1)
        );
        v_G = _v[x_G][y_G];
    }
    else {
        u_G = 0.5 * (_u[x_G][y_G] + _u[x_G+1][y_G]);
        v_G = 0.5 * (_v[x_G][y_G] + _v[x_G][y_G+1]);
    }
}

// TODO: cell label
double SmokeSolver2D::u_with_bnd(int x, int y) {
    if (x >= 0 && x <= _nx && y >= 0 && y < _ny) {
        return _u[x][y];
    } else { // Boundary condition:
        return _wind_u;
    }
}

double SmokeSolver2D::v_with_bnd(int x, int y) {
    if (x >= 0 && x < _nx && y >= 0 && y <= _ny) {
        return _v[x][y];
    } else { // Boundary condition:
        return _wind_v;
    }
}

double SmokeSolver2D::T_with_bnd(int x, int y) {
    if (x >= 0 && x < _nx && y >= 0 && y < _ny) {
        return _T[x][y];
    } else { // Boundary condition:
        return _amb_T;
    }
}

double SmokeSolver2D::s_with_bnd(int x, int y) {
    if (x >= 0 && x < _nx && y >= 0 && y < _ny) {
        return _s[x][y];
    } else { // Boundary condition:
        return _amb_s;
    }
}

CellLabel SmokeSolver2D::label(int x, int y) {
    if (x >= 0 && x < _nx && y >= 0 && y < _ny) {
        return _label[x][y];
    } else { // Boundary condition
        // Assuming surrounded by air
        return FLUID;
    }
}