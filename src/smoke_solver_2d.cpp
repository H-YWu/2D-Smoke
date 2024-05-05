#include "smoke_solver_2d.h"
#include "interpolation.h"
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
    int grid_width, int grid_height, double dx,
    double alpha, double beta,
    double ambient_T, double ambient_s,
    double wind_u, double wind_v,
    double rate_T, double rate_s, double T_target,
    IntegrationScheme integration_scheme
)
    : _nx(grid_width), _ny(grid_height), _dx(dx)
    , _alpha(alpha), _beta(1.0/ambient_T)
    , _amb_T(ambient_T), _amb_s(ambient_s)
    , _wind_u(wind_u), _wind_v(wind_v)
    , _r_T(rate_T), _r_s(rate_s), _T_tar(T_target)
    , _sources(0)
    , _label(grid_width, std::vector<CellLabel>(grid_height))
    , _pressure(grid_width, std::vector<double>(grid_height))
    , _u(grid_width+1, std::vector<double>(grid_height))
    , _v(grid_width, std::vector<double>(grid_height+1))
    , _u_nxt(grid_width+1, std::vector<double>(grid_height))
    , _v_nxt(grid_width, std::vector<double>(grid_height+1))
    , _T(grid_width, std::vector<double>(grid_height))
    , _s(grid_width, std::vector<double>(grid_height))
    , _T_nxt(grid_width, std::vector<double>(grid_height))
    , _s_nxt(grid_width, std::vector<double>(grid_height))
    , _s_max(_amb_s), _dT_max(0.0), _u_max(wind_u), _v_max(wind_v)
    , _intsch(integration_scheme)
{
    init();
}

void SmokeSolver2D::init() {
    // Smoke sources
    int xc = static_cast<int>(0.5 * _nx);
    for (int i = 0; i <= 20; i ++) {
        for (int j = 0; j <= 3; j ++) {
            // TODO: check valid index
            _sources.emplace_back(std::make_pair(xc+i, j));
            if (i != 0)
            _sources.emplace_back(std::make_pair(xc-i, j));
        }
    }

    for (int x = 0; x < _nx; x ++) {
        for (int y = 0; y < _ny; y ++) {
            _label[x][y] = FLUID;
            _T[x][y] = _amb_T; 
            _s[x][y] = _amb_s;
            _T_nxt[x][y] = _amb_T; 
            _s_nxt[x][y] = _amb_s;
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
    double dt = 0.01;
    CFL_dt(dt);

    advect(U, dt);
    advect(V, dt);
    _u.swap(_u_nxt);
    _v.swap(_v_nxt);

    // diffuse();
    force(dt);
    project(dt);

    advect(ST, dt);
    advect(SS, dt);
    _T.swap(_T_nxt);
    _s.swap(_s_nxt);

    step_source(dt);
}

void SmokeSolver2D::CFL_dt(double &dt) {
    int C = 5;
    double dt_max_x = dt, dt_max_y = dt;
    double acc1 = std::abs((_alpha * _s_max) * _g), acc2 = std::abs((_beta * _dT_max) * _g);
    if (_u_max > 0) dt_max_x = std::max(dt_max_x, (C*_dx)/_u_max);
    _v_max = std::max(_v_max+sqrt(C*_dx*acc1), _v_max+sqrt(C*_dx*acc2));
    if (_v_max > 0) dt_max_y = std::max(dt_max_y, (C*_dx)/_v_max);
    dt = std::min(dt_max_x, dt_max_y);
}

void SmokeSolver2D::step_source(double dt) {
    for (long unsigned int i = 0; i < _sources.size(); i ++) {
        int x = _sources[i].first;
        int y = _sources[i].second;
        _T[x][y] = _T[x][y] + (1.0 - std::exp(-_r_T * dt)) * (_T_tar - _T[x][y]);
        _s[x][y] = std::min(_s[x][y] + _r_s * dt, 1.0);
    }
}

void SmokeSolver2D::advect(FieldType ft, double dt) {
    int nx = _nx, ny = _ny;
    if (ft == U) nx ++;
    else if (ft == V) ny ++;

    for (int x = 0; x < nx; x ++) {
        for (int y = 0; y < ny; y ++) {
            double x_P, y_P;
            // Default: Backward Euler
            if (_intsch == BACKWARD_EULER)
                backward_Euler(ft, x, y, dt, x_P, y_P);
            else
                RK2(ft, x, y, dt, x_P, y_P);
            // TODO: interpolation method select
            if (ft == U || ft == V) {
                double u_P, v_P;
                blerp_uv(x_P, y_P, u_P, v_P);
                if (ft == U) _u_nxt[x][y] = u_P;
                else if (ft == V) _v_nxt[x][y] = v_P;
            }
            else {
                int x_high = static_cast<int>(std::ceil(x_P));
                int y_high = static_cast<int>(std::ceil(y_P));
                // Boundary condition
                if (x_P < 0 || x_P > _nx-1 || y_P < 0 || y_P > _ny-1) {
                    if (ft == ST)
                        _T_nxt[x][y] = T_with_bnd(x_high, y_high);
                    else if (ft == SS)
                        _s_nxt[x][y] = s_with_bnd(x_high, y_high);
                }
                else {
                    double tx = (x_P-static_cast<double>(x_high-1))/1.0;
                    double ty = (y_P-static_cast<double>(y_high-1))/1.0;
                    if (ft == ST)
                        _T_nxt[x][y] = blerp(
                            _T[x_high-1][y_high-1], _T[x_high][y_high-1],
                            _T[x_high][y_high-1], _T[x_high][y_high],
                            tx, ty);
                    else if (ft == SS)
                        _s_nxt[x][y] = blerp(
                            _s[x_high-1][y_high-1], _s[x_high][y_high-1],
                            _s[x_high][y_high-1], _s[x_high][y_high],
                            tx, ty);
                }
            }
        }
    }
}

void SmokeSolver2D::force(double dt) {
    _s_max = 0.0, _dT_max = 0.0;
    // Gravity/Buoyancy only applies to v
    for (int x = 0; x < _nx; x ++) {
        for (int y = 0; y <= _ny; y ++) {
            double s = 0.5 * (s_with_bnd(x, y-1) + s_with_bnd(x, y));
            double T = 0.5 * (T_with_bnd(x, y-1) + T_with_bnd(x, y));
            double b = 
                (_alpha * s - _beta * (T - _amb_T)) * _g;
            _v_nxt[x][y] = _v[x][y] + dt * b;

            _s_max = std::max(_s_max, s);
            _dT_max = std::max(_dT_max, (T-_amb_T));
        }
    }
    _v.swap(_v_nxt);
}

void SmokeSolver2D::project(double dt) {
    solve_pressure(dt); 
    // Update velocity of fluid from solved pressure
    update_uv_incompressible(dt);
}

void SmokeSolver2D::backward_Euler(FieldType ft, int x_G, int y_G, double dt, double &x_P, double &y_P) {
    double u_G, v_G;
    cell_uv(ft, x_G, y_G, u_G, v_G);

    x_P = x_G - dt * u_G;
    y_P = y_G - dt * v_G;
}

void SmokeSolver2D::RK2(FieldType ft, int x_G, int y_G, double dt, double &x_P, double &y_P) {
    double u_G, v_G;
    cell_uv(ft, x_G, y_G, u_G, v_G);

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
            int r = x*_ny+y;
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
                    A.add_to_element(r, r-_ny, -1.0 * scale_lhs);
                } // EMPTY just 0
                    // j = (x+1, y)
                if (label(x+1,y)==SOLID) {
                    A.add_to_element(r, r, -1.0 * scale_lhs);
                    // u(x+1,y) - usolid(x+1, y)
                    // Assume solids are all still
                    b[r] += scale_rhs * (_u[x+1][y] - 0.0);
                } 
                else if (label(x+1,y)==FLUID) {
                    A.add_to_element(r, r+_ny, -1.0 * scale_lhs);
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
            _pressure[x][y] = result[x*_ny+y];
        }
    }
    // TODO: check solver
}

void SmokeSolver2D::update_uv_incompressible(double dt) {
    _u_max = 0.0, _v_max = 0.0;

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
                        _u[x][y] -= scale * 
                            (p_with_bnd(x,y) - p_with_bnd(x-1,y));
                    }
                } else {
                    // u(x,y) is unknown
                }

                _u_max = std::max(_u_max, _u[x][y]);
            }

            if (x != _nx) { // update v
                if (label(x,y-1)==FLUID || label(x,y)==FLUID) {
                    if (label(x,y-1)==SOLID || label(x,y)==SOLID) {
                        // v(x,y) = vsolid(x,y);
                        // Assume solids are all still
                        _v[x][y] = 0.0;
                    } else {
                        _v[x][y] -= scale * 
                            (p_with_bnd(x,y) - p_with_bnd(x,y-1));
                    }
                } else {
                    // v(x,y) is unknown
                }

                _v_max = std::max(_v_max, _v[x][y]);
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

void SmokeSolver2D::cell_uv(FieldType ft, int x, int y, double &u, double &v) {
    if (ft == U) {
        u = _u[x][y];
        v = 0.25 * (
            v_with_bnd(x-1, y) + v_with_bnd(x-1, y+1) +
            v_with_bnd(x, y) + v_with_bnd(x, y+1)
        );
    }
    else if (ft == V) {
        u = 0.25 * (
            u_with_bnd(x-1, y) + u_with_bnd(x-1, y+1) +
            u_with_bnd(x, y) + u_with_bnd(x, y+1)
        );
        v = _v[x][y];
    }
    else {
        u = 0.5 * (_u[x][y] + _u[x+1][y]);
        v = 0.5 * (_v[x][y] + _v[x][y+1]);
    }
}

// TODO: cell label SOLID EMPTY extrapolate

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

double SmokeSolver2D::p_with_bnd(int x, int y) {
    if (x >= 0 && x < _nx && y >= 0 && y < _ny) {
        return _pressure[x][y];
    } else { // Boundary condition:
        // Free air pressure
        return 101.325; // Pa
        // return 0.0;
    }
}

CellLabel SmokeSolver2D::label(int x, int y) {
    if (x >= 0 && x < _nx && y >= 0 && y < _ny) {
        return _label[x][y];
    } else { // Boundary condition
        return EMPTY;
    }
}