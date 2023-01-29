#include "AdamsMoulton.hpp"
#include <cmath>
#include <complex>
#include <fstream>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>

// Note: This example uses GSL (GNU Scientific Libraries) for evaluating `exact'
// bessel functions; only required for testing the result.

const std::string description{R"(
...
)"};

//==============================================================================
int main() {

  std::cout << description << "\n";

  double r0 = 5.0e-4;
  double dt = 0.02;
  int n_steps = 500;
  std::vector<double> r;
  std::vector<double> drdt;
  for (int i = 0; i < n_steps; ++i) {
    auto tr = r0 * (std::exp(i * dt));
    r.push_back(tr);
    drdt.push_back(tr); // for this grid, dr/dt = r
  }

  // Define the DerivativeMatrix for the Bessel equation
  struct BesselDerivative
      : AdamsMoulton::DerivativeMatrix<std::size_t, double> {
    double E;
    std::vector<double> r;
    std::vector<double> drdt;
    BesselDerivative(double tE, std::vector<double> tr,
                     std::vector<double> tdrdt)
        : E(tE), r(tr), drdt(tdrdt) {}
    double a(std::size_t) const final { return 0.0; }
    double b(std::size_t i) const final { return 1.0 * drdt[i]; }
    double c(std::size_t i) const final {
      return -2.0 * drdt[i] * (E + 1.0 / r[i]);
    }
    double d(std::size_t) const final { return 0.0; }
  };

  const double En = -0.5;
  BesselDerivative D{En, r, drdt};

  // Construct the solver:
  AdamsMoulton::ODESolver_2x2<12, std::size_t, double> ode{dt, &D};

  double f0 = 2.0 * r0;
  double g0 = (-1.0 + 1 / r0) * f0;
  std::cout << "Initial points: t0=" << r0 << ", f0=" << f0 << ", g0=" << g0
            << "\n";
  std::cout << "Using a K=" << ode.K_steps() << " AM method\n";

  ode.solve_initial_K(0ul, f0, g0);

  // Solve the ODE out to very large t
  std::cout << "...\n";
  std::cout << "Solve ODE all the way out to r=" << r.back()
            << ", with N=" << r.size() << " steps\n";
  std::ofstream out{"out.txt"};
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    out << r.at(i) << " " << ode.f.at(i) << "\n";
  }
  while (ode.last_t() < r.size() - 1) {
    ode.drive();
    out << r.at(ode.last_t()) << " " << ode.last_f() << "\n";
  }
}