#include "euler.h"

// Conservative to primitive
//////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::VectorXd c2p(Eigen::VectorXd Q)
{
    Eigen::VectorXd U(4);
    U(0) = Q(0);
    U(1) = Q(1) / Q(0);
    U(2) = Q(2) / Q(0);
    U(3) = (_gamma - 1.0) * (Q(3) - 0.5 * U(0) * (U(1) * U(1) + U(2) * U(2)));
    return U;
}
// Primitive to conservative
//////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::VectorXd p2c(Eigen::VectorXd U)
{
    Eigen::VectorXd Q(4);
    Q(0) = U(0);
    Q(1) = U(0) * U(1);
    Q(2) = U(0) * U(2);
    Q(3) = U(3) / (_gamma - 1.0) + 0.5 * U(0) * (U(1) * U(1) + U(2) * U(2));
    return Q;
}

//Initial conditions
//////////////////////////////////////////////////////////////////////////////////////////////////////
void initial_conditions(euler &e, Eigen::VectorXd &Q)
{
    Eigen::VectorXd U(4);

    Q = Eigen::VectorXd::Zero(4 * e.N_cells);

    for (int i = 0; i < e.N_cells; i++)
    {
        U(0) = e.rho_inf;
        U(1) = e.u_inf;
        U(2) = e.v_inf;
        U(3) = e.p_inf;

        e.cells[i].Q = p2c(U);

        Q(4 * i) = e.cells[i].Q(0);
        Q(4 * i + 1) = e.cells[i].Q(1);
        Q(4 * i + 2) = e.cells[i].Q(2);
        Q(4 * i + 3) = e.cells[i].Q(3);
    }
}

// Max
//////////////////////////////////////////////////////////////////////////////////////////////////////
double max(double a, double b)
{
    if (a > b)
    {
        return a;
    }
    return b;
}
// Min
//////////////////////////////////////////////////////////////////////////////////////////////////////
double min(double a, double b)
{
    if (a > b)
    {
        return b;
    }
    return a;
}