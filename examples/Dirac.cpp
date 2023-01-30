#include "AdamsMoulton.hpp"
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>

const std::string description{R"(
Solve the radial Dirac equation for Hydrogen (n=1, kappa=-1).
[See W. R. Johnson, Atomic Structure Theory (2007) for details]

In all other examples, the system of ODEs is a single second-order ODE.
The Dirac equation, on the other hand, involves directly solving a pair of 
coupled first-order ODEs.

The Dirac equation can be written as:
    dF/dr = D(t) * F

where:
    F(t) = ( f(r) )
           ( g(r) ),
and (for kappa=-1):
    D(t) = ( -ck/r   (E-V+2c^2)) ) * alpha
           ( -(E-V)  ck/r        )         .

alpha =~ 1/137 is the fine structure constant, c=1/alpha, V = -1/r, and k=-1.

In this example, we take E = +1.0 (i.e., a continuum state).
We don't normalise the function.
)"};

//==============================================================================
int main() {

  std::cout << description << "\n";

  // Set up the non-linear (logarithmically-spaced) radial grid:
  double r0 = 1.0e-3;
  double dr = 0.01;

  // Define the DerivativeMatrix for the Dirac equation
  struct DiracDerivative : AdamsMoulton::DerivativeMatrix<double, double> {
    double E;
    double alpha = 1.0 / 137.036; // approx
    double twoc2 = 2.0 / alpha / alpha;
    DiracDerivative(double tE) : E(tE) {}
    double a(double r) const final { return +1.0 / r; }
    double b(double r) const final { return (E + 1.0 / r + twoc2) * alpha; }
    double c(double r) const final { return -(E + 1.0 / r) * alpha; }
    double d(double r) const final { return -1.0 / r; }
  };

  const double En = +1.0;
  DiracDerivative D{En};

  // Construct the solver:
  AdamsMoulton::ODESolver2D<6, double, double> ode{dr, &D};

  // Set initial points (See Johnson 2007)
  const auto gamma = std::sqrt(1.0 - D.alpha * D.alpha);
  double f0 = 2.0 * std::pow(r0, gamma);
  double g0 = -D.alpha / (gamma + 1) * f0;
  std::cout << "Initial points: t0=" << r0 << ", f0=" << f0 << ", g0=" << g0
            << "\n";

  // Automatically solve the first K points
  ode.solve_initial_K(r0, f0, g0);

  // Solve the ODE out to end of the grid
  double r_max = 50.0;
  std::cout << "...\n";
  std::cout << "Solve ODE all the way out to r=" << r_max << "\n";

  std::string ofname{"Dirac.txt"};
  std::ofstream out{ofname};
  out << "r     f(r)    g(r)\n";
  double r = r0;
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    out << r << " " << ode.f.at(i) << " " << ode.g.at(i) << "\n";
    r += dr;
  }
  while (ode.last_t() < r_max) {
    ode.drive();
    out << ode.last_t() << " " << ode.last_f() << " " << ode.last_g() << "\n";
  }
  std::cout << "Solution written to " << ofname << " for plotting\n";
}