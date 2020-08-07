
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <regex>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

#include "edge.h"
#include "cell.h"

//Solver constants
const double _g = 9.80665;
const double _gamma = 1.4;

//Solver structure
struct euler
{
    std::fstream file;
    std::string case_name;
    double dt;
    int N_iterations;
    int N_cells;
    int N_edges;
    std::vector<cell> cells;
    int output_frequency;
    int riemann_solver;
    int integrator;
    std::vector<edge> edges;
    std::vector<Eigen::Vector3d> vertices;
    int N_vertices;
    int implicit;
    double rho_inf;
    double u_inf;
    double v_inf;
    double p_inf;
    double CFL;
    double M_inf;
    double alpha;
    double tolerance;
};

//Solver functions

//read_config.cpp

void read_config(euler &);

//read_grid.cpp
void read_grid(euler &);

//write_results.cpp
void write_results(euler &, int);

//residual.cpp
Eigen::VectorXd flux(Eigen::VectorXd, Eigen::VectorXd, Eigen::Vector3d, int, double &);

//residual.cpp
Eigen::VectorXd residual(euler &, Eigen::VectorXd &);
Eigen::VectorXd residual_norm(euler &, Eigen::VectorXd &);

//helpers.cpp
Eigen::VectorXd c2p(Eigen::VectorXd);
Eigen::VectorXd p2c(Eigen::VectorXd);
void initial_conditions(euler &, Eigen::VectorXd &);
double max(double, double);
double min(double, double);

//misc.cpp
void local_time_step(euler &);
Eigen::VectorXd rhs(euler &, Eigen::VectorXd);
Eigen::VectorXd freestream(double, double, double, double);
Eigen::VectorXd slip_wall(Eigen::VectorXd Ql, Eigen::Vector3d n);
