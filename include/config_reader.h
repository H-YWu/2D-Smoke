#ifndef CONFIG_READER_H
#define CONFIG_READER_H

#include "yaml-cpp/yaml.h"
#include "smoke_solver_2d.h"

struct SolverConfig {
    int width, height;
    double dx;
    double alpha, beta;
    double amb_T, amb_s, wind_u, wind_v;
    double rate_T, rate_s, T_target;
    IntegrationScheme int_sch;
};


SolverConfig read_YAML(const std::string& filename) {
    SolverConfig config;

    YAML::Node config_file = YAML::LoadFile(filename);

    config.width = config_file["grid"]["width"].as<int>();
    config.height = config_file["grid"]["height"].as<int>();
    config.dx = config_file["grid"]["dx"].as<double>();

    config.alpha = config_file["buoyancy"]["alpha"].as<double>();
    config.beta = config_file["buoyancy"]["beta"].as<double>();

    config.amb_T = config_file["air"]["tempreature"].as<double>();
    config.amb_s = config_file["air"]["concentration"].as<double>();
    config.wind_u = config_file["air"]["u"].as<double>();
    config.wind_v = config_file["air"]["v"].as<double>();

    config.rate_T = config_file["smoke"]["tempreature_rate"].as<double>();
    config.rate_s = config_file["smoke"]["concentration_rate"].as<double>();
    config.T_target = config_file["smoke"]["target_tempreature"].as<double>();

    std::string int_sch_string = config_file["scheme"]["integration"].as<std::string>();
    if (int_sch_string == "BACKWARD_EULER")
        config.int_sch = BACKWARD_EULER;
    else
        config.int_sch = RK2;

    return config;
}

#endif // CONFIG_READER_H