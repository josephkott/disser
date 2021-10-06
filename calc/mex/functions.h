#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <cmath>

#include "Point.h"
using namespace NR;
using namespace std;

// NLS with piecewise constant nonlinearity and absence ordinary (trap) potential.
// Equation for nonlinear stationary modes u(x) e^{i \omega t}:
// $u_{xx} - \omega u + \sigma(x) u^3 = 0$,
// where $\sigma(x) = -1, x \in [0, L_1)$ & $\sigma = +1, x \in [L_1, L_2)$,
// $L = L_1 + L_2$ -- period of nonlinear potential.
// Parameters: $[\omega L_1 L_2]$.
Point<2> piecewise(double x, Point<2> u, double* param)
{
	Point<2> du;
	double omega = param[0], L1 = param[1], L2 = param[2];
    double L = L1 + L2;
    
    double n = floor(x / L);
    double x0 = x - n * L;
    double sigma = (x0 < L1 ? -1 : +1);
    
	du[0] = u[1];
	du[1] = omega * u[0] - sigma * pow(u[0], 3);
	return du;
}

// NLS with cosine pseudopotential and absence of ordinary (trap) potential.
// Equation for nonlinear stationary modes u(x) e^{i \omega t}:
// $u_{xx} - {\omega} u + (\alpha + \cos{2x}) u^3 = 0$.
// Parameters: $[\omega \alpha]$.
Point<2> cosine(double x, Point<2> u, double* param)
{
	Point<2> du;
	double omega = param[0], alpha = param[1];

	du[0] = u[1];
	du[1] = omega * u[0] - (alpha + cos(2 * x)) * pow(u[0], 3);
	return du;
}

// NLS with cosine pseudopotential and harmonic oscillator potential.
// Equation for nonlinaer stationary modes u(x) e^{-i \omega t}:
// $u_{xx} + (\omega - x^2) + (\alpha + \beta \cos Tx) u^3$.
// Parameters: $[\omega T \alpha \beta]$.
Point<2> cosine_nho(double x, Point<2> u, double* param)
{
	Point<2> du;
	double omega = param[0], T = param[1], alpha = param[2], beta = param[3];
	
	du[0] = u[1];
	du[1] = -(omega - pow(x, 2)) * u[0] - (alpha + beta * cos(T * x)) * pow(u[0], 3);
	return du;
}


#endif