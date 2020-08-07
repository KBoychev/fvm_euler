#include "euler.h"

// Flux
Eigen::VectorXd flux(Eigen::VectorXd Ql, Eigen::VectorXd Qr, Eigen::Vector3d n, int riemann_solver, double &lambda_max)
{

	double rhol, ul, vl, El, pl, Hl, Vnl;
	double rhor, ur, vr, Er, pr, Hr, Vnr;
	double nx, ny, nz;
	Eigen::VectorXd deltaQ, Fxl, Fxr, Fyl, Fyr, Fl, Fr, F, w;

	nx = n(0);
	ny = n(1);
	nz = n(2);

	// Difference between right and left conservative variables
	deltaQ = Qr - Ql;

	// Left primitive variables
	rhol = Ql(0);
	ul = Ql(1) / Ql(0);
	vl = Ql(2) / Ql(0);
	El = Ql(3) / Ql(0);
	pl = (_gamma - 1.0) * rhol * (El - 0.5 * (ul * ul + vl * vl));
	Hl = El + pl / rhol;
	Vnl = ul * nx + vl * ny;

	// Right primitive variables
	rhor = Qr(0);
	ur = Qr(1) / Qr(0);
	vr = Qr(2) / Qr(0);
	Er = Qr(3) / Qr(0);
	pr = (_gamma - 1.0) * rhor * (Er - 0.5 * (ur * ur + vr * vr));
	Hr = Er + pr / rhor;
	Vnr = ur * nx + vr * ny;

	//Physical fluxes

	Fxl = Eigen::VectorXd::Zero(4);
	Fxr = Eigen::VectorXd::Zero(4);
	Fyl = Eigen::VectorXd::Zero(4);
	Fyr = Eigen::VectorXd::Zero(4);

	Fxl(0) = rhol * ul;
	Fxl(1) = rhol * ul * ul + pl;
	Fxl(2) = rhol * ul * vl;
	Fxl(3) = rhol * ul * Hl;

	Fxr(0) = rhor * ur;
	Fxr(1) = rhor * ur * ur + pr;
	Fxr(2) = rhor * ur * vr;
	Fxr(3) = rhor * ur * Hr;

	Fyl(0) = rhol * vl;
	Fyl(1) = rhol * vl * ul;
	Fyl(2) = rhol * vl * vl + pl;
	Fyl(3) = rhol * vl * Hl;

	Fyr(0) = rhor * vr;
	Fyr(1) = rhor * vr * ur;
	Fyr(2) = rhor * vr * vr + pr;
	Fyr(3) = rhor * vr * Hr;

	Fl = Fxl * nx + Fyl * ny;
	Fr = Fxr * nx + Fyr * ny;

	// Lax-Friedriechs flux
	//////////////////////////////////////////////////////////////////////

	if (riemann_solver == 0)
	{
		double al, ar, lambdal_max, lambdar_max;
		al = std::sqrt(_gamma * pl / rhol);
		ar = std::sqrt(_gamma * pr / rhor);
		lambdal_max = std::abs(Vnl) + al;
		lambdar_max = std::abs(Vnr) + ar;
		lambda_max = max(lambdal_max, lambdar_max);
		w = lambda_max * deltaQ;
	}

	//Roe flux
	//////////////////////////////////////////////////////////////////////

	if (riemann_solver == 1)
	{
		double C, rho, u, v, H, a, Vn;
		double drho, du, dv, dp, dVn;
		Eigen::VectorXd I(4), lambda(4), r1(4), r2(4), r3(4), r4(4);

		C = std::sqrt(rhor / rhol);
		rho = C * rhol;
		u = (ul + C * ur) / (1.0 + C);
		v = (vl + C * vr) / (1.0 + C);
		H = (Hl + C * Hr) / (1.0 + C);
		a = std::sqrt((_gamma - 1.0) * (H - 0.5 * (u * u + v * v)));
		Vn = u * nx + v * ny;

		drho = rhor - rhol;
		du = ur - ul;
		dv = vr - vl;
		dp = pr - pl;
		dVn = Vnr - Vnl;

		I(0) = (dp - rho * a * dVn) / (2.0 * a * a);
		I(1) = (dp + rho * a * dVn) / (2.0 * a * a);
		I(2) = drho - dp / (a * a); //here is the problem
		I(3) = rho;

		lambda(0) = std::abs(Vn - a);
		lambda(1) = std::abs(Vn + a);
		lambda(2) = std::abs(Vn);
		lambda(3) = std::abs(Vn);

		for (int i = 0; i < 4; i++)
		{
			double dlambda;
			dlambda = 0.2 * a;
			if (lambda(i) < dlambda)
			{
				lambda(i) = 0.5 * (lambda(i) * lambda(i) / dlambda + dlambda);
			}
		}

		r1(0) = 1.0;
		r1(1) = u - a * nx;
		r1(2) = v - a * ny;
		r1(3) = H - a * Vn;

		r2(0) = 1.0;
		r2(1) = u + a * nx;
		r2(2) = v + a * ny;
		r2(3) = H + a * Vn;

		r3(0) = 1.0;
		r3(1) = u;
		r3(2) = v;
		r3(3) = 0.5 * (u * u + v * v);

		r4(0) = 0.0;
		r4(1) = du - dVn * nx;
		r4(2) = dv - dVn * ny;
		r4(3) = u * du + v * dv - Vn * dVn;

		w = lambda(0) * I(0) * r1 + lambda(1) * I(1) * r2 + lambda(2) * I(2) * r3 + lambda(3) * I(3) * r4;

		lambda_max = std::abs(Vn) + a;
	}

	F = 0.5 * (Fl + Fr) - 0.5 * w;

	return F;
}
