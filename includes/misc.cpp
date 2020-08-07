#include "euler.h"

//Local time step
//////////////////////////////////////////////////////////////////////////////////////////////////////

void local_time_step(euler &e)
{
    for (int i = 0; i < e.N_cells; i++)
    {
        e.cells[i].dtau = e.CFL * e.cells[i].S / (0.5 * e.cells[i].lambda_max);
    }
}

// Right hand side
//////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::VectorXd rhs(euler &e, Eigen::VectorXd R)
{
    Eigen::VectorXd rhs;
    rhs = Eigen::VectorXd::Zero(R.rows());

    //Set residual vector
    for (int i = 0; i < e.N_cells; i++)
    {
        rhs(4 * i) = -e.cells[i].dtau / e.cells[i].S * R(4 * i);
        rhs(4 * i + 1) = -e.cells[i].dtau / e.cells[i].S * R(4 * i + 1);
        rhs(4 * i + 2) = -e.cells[i].dtau / e.cells[i].S * R(4 * i + 2);
        rhs(4 * i + 3) = -e.cells[i].dtau / e.cells[i].S * R(4 * i + 3);
    }

    return rhs;
}

// Freestream boundary condition
//////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::VectorXd freestream(double rho_inf, double u_inf, double v_inf, double p_inf)
{
    Eigen::VectorXd Ub(4);
    Ub(0) = rho_inf;
    Ub(1) = u_inf;
    Ub(2) = v_inf;
    Ub(3) = p_inf;
    return p2c(Ub);
}

// Slip wall boundary condition
//////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::VectorXd slip_wall(Eigen::VectorXd Ql, Eigen::Vector3d n)
{
    double rhol, ul, vl, El;
    double nx, ny, nz;
    double Vn;
    Eigen::VectorXd Ul, Ub;

    Ul = c2p(Ql);
    nx = n(0);
    ny = n(1);

    Vn = Ul(1) * nx + Ul(2) * ny;

    Ub = Ul;

    Ub(1) = Ub(1) - Vn * nx;
    Ub(2) = Ub(2) - Vn * ny;

    return p2c(Ub);
}