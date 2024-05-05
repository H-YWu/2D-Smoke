#ifndef CONFIG_H
#define CONFIG_H

// Default parameters for smoke simulation
int grid_w = 1280;
int grid_h = 720;
double dx = 1.0;
double alpha=0.5, beta=0.1;
double ambient_T=273.0, ambient_s=0.0;
double wind_u=1.0, wind_v=0.0;
double rate_T=8.0, rate_s=10.0, T_target=300.0;

#endif // CONFIG_H 