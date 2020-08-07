#include "euler.h"

//Residual  norm
//////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::VectorXd residual_norm(euler &e, Eigen::VectorXd &R)
{
    Eigen::VectorXd R_norm;
    R_norm = Eigen::VectorXd::Zero(4);
    for (int i = 0; i < e.N_cells; i++)
    {
        R_norm(0) += std::abs(R(4 * i + 0));
        R_norm(1) += std::abs(R(4 * i + 1));
        R_norm(2) += std::abs(R(4 * i + 2));
        R_norm(3) += std::abs(R(4 * i + 3));

        if (std::isnan(R_norm(0)) || std::isnan(R_norm(1)) || std::isnan(R_norm(2)) || std::isnan(R_norm(3)))
        {
            std::cout << "Residual in cell " << i << " is nan, stopping..." << std::endl;
            std::exit(-1);
        }
    }
    R_norm /= e.N_cells;
    return R_norm;
}

// Residual
//////////////////////////////////////////////////////////////////////////////////////////////////////
//
//          Edge
//            o
//            | Normal
// Left cell  |------> n Right cell
//    cl      |              cr
//            o
//
// 1. Compute the numerical flux
// 2. Subtract it to the residual of cell 1 and add it to the residual of c2
//
// Note: If the edge is at a boundary there is no right cell!
//
//////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::VectorXd residual(euler &e, Eigen::VectorXd &Q)
{
    Eigen::VectorXd R;
    R = Eigen::VectorXd::Zero(Q.rows());

    // Set cell conservative variables and reset cell residuals
    for (int i = 0; i < e.N_cells; i++)
    {
        e.cells[i].Q(0) = Q(4 * i);
        e.cells[i].Q(1) = Q(4 * i + 1);
        e.cells[i].Q(2) = Q(4 * i + 2);
        e.cells[i].Q(3) = Q(4 * i + 3);
        e.cells[i].R(0) = 0;
        e.cells[i].R(1) = 0;
        e.cells[i].R(2) = 0;
        e.cells[i].R(3) = 0;
        e.cells[i].lambda_max = 0;
    }

    // Set left and right conservative variables, apply boundary conditions
    for (int i = 0; i < e.N_edges; i++)
    {
        int cl, cr;
        Eigen::Vector3d n;
        Eigen::VectorXd Ql, Qr;

        cl = e.edges[i].celll; //Left cell of edge
        cr = e.edges[i].cellr; //Right cell of edge
        n = e.edges[i].n;      //Edge normal

        Ql = e.cells[cl].Q; // Left conservative variables

        if (cr != -1)
        {
            Qr = e.cells[cr].Q; //Right conservative variables if right cell exists
        }
        else
        {
            //Booundary condition if right cell does not exist

            if (e.cells[cl].type == 1) //Freestream boundary condition
            {
                Qr = freestream(e.rho_inf, e.u_inf, e.v_inf, e.p_inf);
            }

            if (e.cells[cl].type == 2) //Slip wall boundary condition
            {
                Qr = slip_wall(Ql, n);
            }
        }

        e.edges[i].Ql = Ql;
        e.edges[i].Qr = Qr;
    }

    // Residual calculation
    for (int i = 0; i < e.N_edges; i++)
    {
        int cl, cr;
        double l, lambda_max;
        Eigen::Vector3d n;
        Eigen::VectorXd Ql, Qr, F;

        cl = e.edges[i].celll; //Left cell of edge
        cr = e.edges[i].cellr; //Right cell of edge
        n = e.edges[i].n;      //Edge normal
        l = e.edges[i].l;      //Edge length

        Ql = e.edges[i].Ql; // Left conservative variables
        Qr = e.edges[i].Qr; // Right conservative variables

        F = flux(Ql, Qr, n, e.riemann_solver, lambda_max); // Flux

        // Add flux to left cell
        e.cells[cl].R += F * l;
        e.cells[cl].lambda_max += lambda_max * l;

        if (cr != -1)
        {
            // Subtract flux from right cell if right cell exists
            e.cells[cr].R -= F * l;
            e.cells[cr].lambda_max += lambda_max * l;
        }
    }

    //Set residual vector
    for (int i = 0; i < e.N_cells; i++)
    {
        R(4 * i) = e.cells[i].R(0);
        R(4 * i + 1) = e.cells[i].R(1);
        R(4 * i + 2) = e.cells[i].R(2);
        R(4 * i + 3) = e.cells[i].R(3);
    }

    return R;
}