# Finite volume method for euler equations

A 2D unstructured finite volume method (FVM) euler solver written in C++. Fluxes can be evaluated with the Lax–Friedrichs or the Roe method. Explicit pseudo-time stepping is available. Expicit pseudo-time stepping can be performed with the Euler, the Runge-Kutta 2nd-order, and the Runge-Kutta 4th-order methods. The Eigen C++ template library for linear algebra is used along with OpenMP. 

Effect of the integration and flux method; Roe and Runge-Kutta 2nd-order (left), Roe and Runge-Kutta 4th-order (middle), Lax-Friedrichs + Euler 1st order (right) 

![Effect of integration and flux method](https://github.com/KBoychev/fvm_euler/blob/master/naca4412_pressure.png "Integration and flux method")

# Compile

The code requires the Eigen C++ template libary. Download the libary and extract its contents to a folder named Eigen inside the includes folder (includes/Eigen). Once extracted use make to compile the code. The sw_solver binary will be created in the bin folder. 

# Run

To run the code navigate to the folder of the case e.g. naca4412 and run ../bin/euler_solver naca4412. The program will solve the flow around a NACA4412 airfoil using the Roe + Euler method, until the residual ratio is less than 1e-5. The full list of options is:


* dt - timestep
* N_iterations - number of iterations
* riemann_solver - Riemman solver; 1 - Lax-Friedrichs, 2 - Roe
* output_frequency - Results output frequency; 0 - no output, i where i > 0 - every ith iteration 
* integrator - Integration method; 0 - Euler, 1 - Runge-Kutta 2nd-order, 2 - Runge-Kutta 4th order
* implicit - Solve the equation explicitly or implicitly; 0 - Explicit, 1 - Implicit (ONLY EXPLICIT IS SUPPORTED AT THE MOMENT)

Note: If the implicit option is set to 1 the integrator option is not used. 

The grid file (.grd) has the following format

```
N_vertices N_cell
x1 y1 z1
x2 y2 z2
xN yN zN
f11 f12 f13
f21 f22 f23
fN1 fN2 fN3
```

Note: Quadrilateral cells can be defined as follows: f11 f12 f13 f14 where f14 is the id of the 4th vertex of the 1st cell.