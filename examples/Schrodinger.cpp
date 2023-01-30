#include "AdamsMoulton.hpp"
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

const std::string description{R"(
Solve the radial Schrodinger equation for Hydrogen (n=1, l=0).
[See W. R. Johnson, Atomic Structure Theory (2007) for details]

Psi_nlm(r,theta,phi) = (1/r) P_nl(r) * Y_lm(theta,phi)

P''(r) + 2[E - V(r) - l(l+1)/(2r^2)]P(r) = 0

We will take v(r) = -1/r, l=0, and E = -0.5 (the lowest bound state).

We use this as an example of a non-linear grid, and the case where the function values for the derivative matrix are stored in arrays*.

We define:

r(t)  = r0 * exp(t)
dr/dt = r0 * exp(t)

and store r[i] on a grid, with 
t_i = i * dt     i=(0,1,..., N-1)
dt = 0.02

The Schrodinger equation can be written as:
    dF/dt = dF/dr * dr/dt = D(t) * F

where:
    F(t) = F[i] = F(r[i]) = ( P(r[i]) )
                            ( dP/dr   )
for the ith point along the grid.

The Jacobian, drdt, must be included in D:

and:
    D(t) = ( 0                      dr/dt[i] )
           ( -2*dr/dt[i]*(E-v[i])   0        ).

*In this particular case, we have a functional form for r[t] and v(r) - but
in more complicated cases we don't, which is why it can be useful to store
them on a grid.

)"};

//==============================================================================
int main() {

  std::cout << description << "\n";

  // Set up the non-linear (logarithmically-spaced) radial grid:
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

  // Define the DerivativeMatrix for the Schrodinger equation
  // Note: it stores the r values as arrays (std::vector), so takes an array
  // index (std::size_t) as the function argument (and returns double).
  // In this particular case, we have a functional form for r[t] and v(r) - but
  // in more complicated cases we don't, which is why it can be useful to store
  // on a grid.
  struct SchrodingerDerivative
      : AdamsMoulton::DerivativeMatrix<std::size_t, double> {
    double E;
    std::vector<double> r;
    std::vector<double> drdt;
    SchrodingerDerivative(double tE, std::vector<double> tr,
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
  SchrodingerDerivative D{En, r, drdt};

  // Construct the solver:
  AdamsMoulton::ODESolver_2x2<12, std::size_t, double> ode{dt, &D};

  // Set initial points (See Johnson 2007)
  double f0 = 2.0 * r0;
  double g0 = (-1.0 + 1 / r0) * f0;
  std::cout << "Initial points: t0=" << r0 << ", f0=" << f0 << ", g0=" << g0
            << "\n";

  // Automatically solve the first K points
  // Note: the initial 't' here is zero, since in this case it is an array index
  // - corresponding to r[0]
  ode.solve_initial_K(0ul, f0, g0);

  // Solve the ODE out to end of the grid
  std::cout << "...\n";
  std::cout << "Solve ODE all the way out to r=" << r.back()
            << ", with N=" << r.size() << " steps\n";

  std::string ofname{"Schrodinger.txt"};
  std::ofstream out{ofname};
  out << "i     r(i)     P(r)\n";
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    out << i << " " << r.at(i) << " " << ode.f.at(i) << "\n";
  }
  while (ode.last_t() < r.size() - 1) {
    ode.drive();
    // last_t() here return the array index, i
    out << ode.last_t() << " " << r.at(ode.last_t()) << " " << ode.last_f()
        << "\n";
  }
  std::cout << "Solution written to " << ofname << " for plotting\n";
}