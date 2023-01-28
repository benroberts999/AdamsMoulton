#include "AdamsMoulton.hpp"
#include <cmath>
#include <complex>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>

//==============================================================================
int main() {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Adams Moulton example: Bessel equation\n";

  const int n = 0;

  // This uses GSL library, only to compare against exact solution
  auto y = [=](double x) { return gsl_sf_bessel_Jn(n, double(x)); };
  auto dy = [=](double x) {
    return 0.5 * (gsl_sf_bessel_Jn(n - 1, double(x)) -
                  gsl_sf_bessel_Jn(n + 1, double(x)));
  };

  // Define the DerivativeMatrix
  struct BesselDerivative : AdamsMoulton::DerivativeMatrix<double, double> {
    int n;
    BesselDerivative(int tn) : n(tn) {}
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double t) const final { return std::pow(n / t, 2) - 1.0; }
    double d(double t) const final { return -1.0 / t; }
  };

  BesselDerivative D{n};

  double dt = 0.001;

  AdamsMoulton::ODESolver_2x2<12, double, double> ode{dt, &D};

  double t0 = 1.0e-4;
  // Note: These are approximate, since f(0)=1, but t0=1.0e-6
  double f0 = 1.0;
  double g0 = 0.0; // dy(t0);
  std::cout << "Initial points: t0=" << t0 << ", f0=" << f0 << ", g0=" << g0
            << "\n";

  ode.solve_initial_K(t0, f0, g0);

  // Print initial points, compare to exact
  double t = t0;
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    printf("%.4f %.10f [%.10f] %.10f [%.10f]\n", t, ode.f[i], y(t), ode.g[i],
           dy(t));
    t += dt;
  }

  // Solve the ODE out to very large t
  std::cout << "...\n";
  while (ode.last_t() < 100.0) {
    ode.drive();
  }

  const auto final_f = ode.last_f();
  const auto final_t = ode.last_t();
  const auto expected_f = y(final_t);

  auto eps = (final_f - expected_f) / expected_f;
  printf("%.4f %.10f [%.10f] %.10f [%.10f]  %.2e\n", ode.last_t(), ode.last_f(),
         y(final_t), ode.last_g(), dy(final_t), eps);
}