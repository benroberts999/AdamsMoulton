#include "AdamsMoulton.hpp"
#include <cmath>
#include <complex>
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

The exact solution, with y(0)=1 and y'(0)=0, is the Bessel function: 
y(t) = J_n(t)
)"};

//==============================================================================
int main() {

  std::cout << description << "\n";

  // Define the DerivativeMatrix for the Bessel equation
  struct BesselDerivative : AdamsMoulton::DerivativeMatrix<double> {
    int n;
    BesselDerivative(int tn) : n(tn) {}
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double t) const final { return std::pow(n / t, 2) - 1.0; }
    double d(double t) const final { return -1.0 / t; }
  };

  const int n = 0;
  BesselDerivative D{n};

  // Exact solution: just to compare result
  // This uses GSL library, only to compare against exact solution
  auto y = [=](double t) { return gsl_sf_bessel_Jn(n, t); };
  auto dy = [=](double t) {
    return 0.5 * (gsl_sf_bessel_Jn(n - 1, t) - gsl_sf_bessel_Jn(n + 1, t));
  };

  // Set the stp size
  double dt = 0.001;
  // Construct the solver:
  AdamsMoulton::ODESolver2D<6> ode{dt, &D};

  // Set the initial conditions. Note: derivative matrix has 1/t, so we cannot
  // start from 0. Instead start from 'small' value:
  double t0 = 1.0e-4;
  double f0 = y(t0);
  double g0 = dy(t0);
  // Note: The actual initial value, y(0)=1, is impossible due to 1/t
  // In realistic scenario, we might not know y(t0) exactly.
  // The following approximation works very well
  // double f0 = 1.0;
  // double g0 = 0.0;
  std::cout << "Initial points: t0=" << t0 << ", f0=" << f0 << ", g0=" << g0
            << "\n";
  std::cout << "Using a K=" << ode.K_steps() << " AM method\n";

  ode.solve_initial_K(t0, f0, g0);

  // Print initial points, compare to exact

  std::cout << "Compare initial K points to expected (exact) solution:\n";
  std::cout << "  t       y(t)         [Exact Soln   ]  dy/dt(t)     [Exact "
               "Soln   ]\n";
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    const auto t = ode.t[i];
    printf("%8.4f %13.10f [%13.10f] %13.10f [%13.10f]\n", t, ode.f[i], y(t),
           ode.g[i], dy(t));
  }

  // Solve the ODE out to very large t
  std::cout << "...\n";
  double t_target = 100.0;
  std::cout << "Solve ODE all the way out to t=" << t_target
            << " and compare:\n";
  while (ode.last_t() < t_target) {
    ode.drive();
  }

  const auto final_f = ode.last_f();
  const auto final_t = ode.last_t();
  const auto expected_f = y(final_t);

  auto eps = (final_f - expected_f) / expected_f;
  printf("%8.4f %13.10f [%13.10f] %13.10f [%13.10f]  error: %.1e\n",
         ode.last_t(), ode.last_f(), y(final_t), ode.last_g(), dy(final_t),
         eps);

  // also: backwards:
  // Construct the solver:
  std::cout << "\nWe can also drive the ODE 'backwards' by setting -ve dt\n";
  AdamsMoulton::ODESolver2D<12> ode_back{-dt, &D};
  ode_back.solve_initial_K(-t0, y(-t0), dy(-t0));
  double backwards_target = -75.0;
  std::cout << "Start close to zero; drive backwards to t=" << backwards_target
            << ":\n";

  // The "previous" K points are stored in order of the solution, regardless of
  // the direction we drive the ODE:
  for (std::size_t i = 0; i < ode_back.f.size(); ++i) {
    const auto t = ode_back.t[i];
    printf("%8.4f %13.10f [%13.10f] %13.10f [%13.10f]\n", t, ode_back.f[i],
           y(t), ode_back.g[i], dy(t));
  }

  while (ode_back.last_t() > backwards_target) {
    ode_back.drive();
  }
  const auto final_f2 = ode_back.last_f();
  const auto final_t2 = ode_back.last_t();
  const auto expected_f2 = y(final_t2);
  auto eps2 = (final_f2 - expected_f2) / expected_f2;
  std::cout << "...\n";
  printf("%8.4f %13.10f [%13.10f] %13.10f [%13.10f]  error: %.1e\n",
         ode_back.last_t(), ode_back.last_f(), y(final_t2), ode_back.last_g(),
         dy(final_t2), eps2);
}