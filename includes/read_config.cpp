#include "euler.h"

//Read config
void read_config(euler &e)
{

    e.file.open(e.case_name + ".conf", std::fstream::in);

    if (e.file.is_open())
    {
        std::string file_line;
        std::smatch file_line_match;
        std::regex file_line_regex;
        file_line_match.empty();

        e.dt = 1.0;
        e.N_iterations = 1;
        e.N_cells = 0;
        e.output_frequency = 0;
        e.riemann_solver = 0;
        e.integrator = 0;
        e.implicit = 0;
        e.rho_inf = 0;
        e.u_inf = 0;
        e.v_inf = 0;
        e.p_inf = 0;
        e.CFL = 0.5;
        e.M_inf = 0;
        e.alpha = 0;

        while (std::getline(e.file, file_line))
        {

            file_line_regex.assign("(.+)(?: = )(.+)");

            if (std::regex_search(file_line, file_line_match, file_line_regex))
            {

                if (file_line_match[1].compare("dt") == 0)
                {
                    e.dt = std::stod(file_line_match[2]);
                    if (e.dt <= 0)
                    {
                        std::cout << "Timestep (dt) must be greater than 0!" << std::endl;
                        std::exit(-1);
                    }
                }
                if (file_line_match[1].compare("N_iterations") == 0)
                {
                    e.N_iterations = std::stoi(file_line_match[2]);
                    if (e.N_iterations < 0)
                    {
                        std::cout << "Number of iterations (N_iterations) must be greater than 0!" << std::endl;
                        std::exit(-1);
                    }
                }
                if (file_line_match[1].compare("riemann_solver") == 0)
                {
                    e.riemann_solver = std::stoi(file_line_match[2]);

                    if (e.riemann_solver < 0 && e.riemann_solver > 1)
                    {
                        std::cout << "Riemann solver (riemann_solver) must be 0 or 1!" << std::endl;
                        std::exit(-1);
                    }
                }
                if (file_line_match[1].compare("output_frequency") == 0)
                {
                    e.output_frequency = std::stoi(file_line_match[2]);

                    if (e.output_frequency < 0)
                    {
                        std::cout << "Output frequency (output_frequency) must be 0 or greater!" << std::endl;
                        std::exit(-1);
                    }
                }
                if (file_line_match[1].compare("integrator") == 0)
                {
                    e.integrator = std::stoi(file_line_match[2]);

                    if (e.integrator != 0 && e.integrator != 1 && e.integrator != 2)
                    {
                        std::cout << "Integrator (integrator) must be 0,1, or 2!" << std::endl;
                        std::exit(-1);
                    }
                }
                if (file_line_match[1].compare("implicit") == 0)
                {
                    e.implicit = std::stoi(file_line_match[2]);

                    if (e.implicit < 0 && e.implicit > 1)
                    {
                        std::cout << "Implicit (implicit) must be 0 or 1!" << std::endl;
                        std::exit(-1);
                    }
                }
                if (file_line_match[1].compare("rho_inf") == 0)
                {
                    e.rho_inf = std::stod(file_line_match[2]);

                    if (e.rho_inf < 0)
                    {
                        std::cout << "Freestream density (rho_inf) must be greater than 0!" << std::endl;
                        std::exit(-1);
                    }
                }
                if (file_line_match[1].compare("u_inf") == 0)
                {
                    e.u_inf = std::stod(file_line_match[2]);
                }
                if (file_line_match[1].compare("v_inf") == 0)
                {
                    e.v_inf = std::stod(file_line_match[2]);
                }
                if (file_line_match[1].compare("p_inf") == 0)
                {
                    e.p_inf = std::stod(file_line_match[2]);

                    if (e.p_inf < 0)
                    {
                        std::cout << "Freestream pressure (p_inf) must be greater than 0!" << std::endl;
                        std::exit(-1);
                    }
                }
                if (file_line_match[1].compare("CFL") == 0)
                {
                    e.CFL = std::stod(file_line_match[2]);

                    if (e.CFL < 0)
                    {
                        std::cout << "CFL (CFL) must be greater than 0!" << std::endl;
                        std::exit(-1);
                    }
                }
                if (file_line_match[1].compare("M_inf") == 0)
                {
                    e.M_inf = std::stod(file_line_match[2]);

                    if (e.M_inf <= 0)
                    {
                        std::cout << "Freestream Mach number (M_inf) must be greater than 0!" << std::endl;
                        std::exit(-1);
                    }
                }
                if (file_line_match[1].compare("alpha") == 0)
                {
                    e.alpha = std::stod(file_line_match[2]);
                }
                if (file_line_match[1].compare("tolerance") == 0)
                {
                    e.tolerance = std::stod(file_line_match[2]);

                    if (e.tolerance < 0)
                    {
                        std::cout << "Tolerance (tolerance) must be greater than 0!" << std::endl;
                        std::exit(-1);
                    }
                }
            }
        }

        e.rho_inf = 1.0;                                      //rho/rho_ref
        e.u_inf = 1.0 * std::cos(e.alpha * EIGEN_PI / 180.0); //u/V_ref
        e.v_inf = 1.0 * std::sin(e.alpha * EIGEN_PI / 180.0); //v/V_ref
        e.p_inf = 1.0 / (_gamma * e.M_inf * e.M_inf);         //p/(gamma*p_ref*M_ref^2)

        e.file.close();
    }
    else
    {
        std::cout << "Failed reading config file!" << std::endl;
        std::exit(-1);
    }
}