#include <iostream>
#include <string>
#include <chrono>

#include "euler.h"

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		std::cout << "Usage: euler_solver <case>" << std::endl;
		std::exit(-1);
	}

	// Variables
	////////////////////////////////////////////////////////////////////

	euler e;											 //euler solver structure
	Eigen::VectorXd Q, Q_tmp, k1, k2, k3, k4, R, R_norm; //solution and residual vectors
	double t, R_norm_max;
	int iteration;

	e.case_name = argv[1];
	t = 0.0;
	iteration = 0;

	// Config
	////////////////////////////////////////////////////////////////////

	std::cout << "Reading config...";

	read_config(e);

	std::cout << "Done!" << std::endl;

	std::cout << std::endl;
	std::cout << "dt = " << e.dt << std::endl;
	std::cout << "N_iterations = " << e.N_iterations << std::endl;
	std::cout << "output_frequency = " << e.output_frequency << std::endl;
	std::cout << "riemann_solver = " << e.riemann_solver << std::endl;
	std::cout << "integrator = " << e.integrator << std::endl;
	std::cout << "implicit = " << e.implicit << std::endl;
	std::cout << "rho_inf = " << e.rho_inf << std::endl;
	std::cout << "u_inf = " << e.u_inf << std::endl;
	std::cout << "v_inf = " << e.v_inf << std::endl;
	std::cout << "p_inf = " << e.p_inf << std::endl;
	std::cout << "M_inf = " << e.M_inf << std::endl;
	std::cout << "alpha = " << e.alpha << std::endl;
	std::cout << "CFL = " << e.CFL << std::endl;
	std::cout << std::endl;

	// Grid
	////////////////////////////////////////////////////////////////////

	std::cout << "Reading grid...";

	read_grid(e);

	std::cout << "Done!" << std::endl;

	std::cout << "N_vertices=" << e.N_vertices << std::endl;
	std::cout << "N_cells=" << e.N_cells << std::endl;
	std::cout << "N_edges=" << e.N_edges << std::endl;

	//Initial conditions
	////////////////////////////////////////////////////////////////////

	std::cout << "Setting initial conditions...";

	initial_conditions(e, Q);

	std::cout << "Done!" << std::endl;

	// Time loop
	////////////////////////////////////////////////////////////////////

	std::cout << "Performing " << e.N_iterations << " iterations..." << std::endl;
	std::cout << "Iteration	max(R)	max(R)/max(R0)" << std::endl;
	std::cout << std::setprecision(3);

	while (true)
	{

		//Residual
		R = residual(e, Q);

		//Residual norm
		R_norm = residual_norm(e, R);

		if (iteration == 0)
		{
			R_norm_max = R_norm.maxCoeff();
		}

		std::cout << std::setprecision(3);
		std::cout << iteration << "	" << R_norm.maxCoeff() << "	" << R_norm.maxCoeff() / R_norm_max << std::endl;

		if (R_norm.maxCoeff() / R_norm_max < e.tolerance)
		{
			write_results(e, iteration);
			std::cout << "Done! Tolerance reached." << std::endl;
			break;
		}

		if (iteration == e.N_iterations)
		{
			write_results(e, iteration);
			std::cout << "Done! Maximum iterations reached." << std::endl;
			break;
		}

		iteration += 1;

		//Local time step
		local_time_step(e);

		//Explicit integration
		if (e.integrator == 0) //Euler 1st-order integration
		{
			Q += rhs(e, R);
		}
		if (e.integrator == 1) //Runge-kutta 2nd-order integration
		{
			k1 = rhs(e, R);
			Q_tmp = Q + k1 / 2.0;
			k2 = rhs(e, residual(e, Q_tmp));
			Q += k2;
		}
		if (e.integrator == 2) //Runge-kutta 4th-order integration
		{

			k1 = rhs(e, R);
			Q_tmp = Q + k1 / 2.0;
			k2 = rhs(e, residual(e, Q_tmp));
			Q_tmp = Q + k2 / 2.0;
			k3 = rhs(e, residual(e, Q_tmp));
			Q_tmp = Q + k3;
			k4 = rhs(e, residual(e, Q_tmp));

			Q += 1.0 / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
		}

		// Results output
		if (e.output_frequency > 0)
		{
			if (iteration % e.output_frequency == 0)
			{
				write_results(e, iteration);
			}
		}
	}

	return 0;
}
