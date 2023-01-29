#include "AdamsMoulton.hpp"
#include <cmath>
#include <complex>
#include <fstream>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>

// Note: This example uses GSL (GNU Scientific Libraries) for evaluating `exact'
// bessel functions; only required for testing the result.

const std::string description{R"(
Adams Moulton example: Bessel equation.

t^2 d^2y/dt^2 + t dy/dt + (t^2-n^2)y = 0

Can be written as:
    dF/dt = D * F

where:
    F(t) = ( y(t)  )
           ( dy/dt )
and:
    D(t) = ( 0                 1 )
           ( (n/t)^2 - 1    -1/t ).

The exact solution, with y(0)=1 and y'(0)=0, is theBessel function: 
y(t) = J_n(t)
)"};

//==============================================================================
int main() {

  std::cout << description << "\n";

  // Define the DerivativeMatrix for the Bessel equation
  struct BesselDerivative : AdamsMoulton::DerivativeMatrix<double> {
    double E;
    BesselDerivative(double tE) : E(tE) {}
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double r) const final { return -2.0 * (E + 1.0 / r); }
    double d(double) const final { return 0.0; }
  };

  const double En = -0.5;
  BesselDerivative D{En};

  // Set the stp size
  double dt = 0.01;
  // Construct the solver:
  AdamsMoulton::ODESolver_2x2<12> ode{dt, &D};

  // Set the initial conditions. Note: derivative matrix has 1/t, so we cannot
  // start from 0. Instead start from 'small' value:
  double t0 = 1.0e-4;
  double f0 = 2.0 * t0;
  double g0 = (-1.0 + 1 / t0) * f0;
  // Note: The actual initial value, y(0)=1, is impossible due to 1/t
  // In realistic scenario, we might not know y(t0) exactly.
  // The following approximation works very well
  // double f0 = 1.0;
  // double g0 = 0.0;
  std::cout << "Initial points: t0=" << t0 << ", f0=" << f0 << ", g0=" << g0
            << "\n";
  std::cout << "Using a K=" << ode.K_steps() << " AM method\n";

  ode.solve_initial_K(t0, f0, g0);

  // Solve the ODE out to very large t
  std::cout << "...\n";
  double t_target = 15.0;
  std::cout << "Solve ODE all the way out to t=" << t_target
            << " and compare:\n";
  auto r = t0;
  std::ofstream out{"out.txt"};
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    out << r << " " << ode.f.at(i) << "\n";
    r += ode.dt();
  }
  while (ode.last_t() < t_target) {
    ode.drive();
    out << ode.last_t() << " " << ode.last_f() << "\n";
  }
}