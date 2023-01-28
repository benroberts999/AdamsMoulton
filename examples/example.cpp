#include "AdamsMoulton.hpp"
#include <cmath>
#include <complex>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>

//==============================================================================
int main() {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Adams Moulton test/example\n";

  /*
  1. Bessel example
    a. few different K
    b. Drive forwards + backwards
    c. double and float ?
    d. Non-uniform grid
  2. Complex example
  3. Array index example
  4. Documentation
  5. Unit tests (innards + outards)
     + CI?
  --
  A. Schro/Dirac example
  B. Add simple ODE case


  */

  int n = 0;

  auto y = [=](double x) { return gsl_sf_bessel_Jn(n, double(x)); };

  auto dy = [=](double x) {
    return 0.5 * (gsl_sf_bessel_Jn(n - 1, double(x)) -
                  gsl_sf_bessel_Jn(n + 1, double(x)));
  };

  struct BesselDerivative : AdamsMoulton::DerivativeMatrix<double, double> {
    int n;
    BesselDerivative(int tn) : n(tn) {}
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double) const final {
      // return std::pow(double{1.0f * n} / t, 2) - 1.0f;
      return -1.0f;
    }
    double d(double t) const final { return -1.0 / t; }
  };

  BesselDerivative D2{n};

  // double dt = {-0.0001, 0.0001};
  double dt = -0.0001f;

  AdamsMoulton::ODESolver_2x2<12, double, double> ode{dt, &D2};

  // double t0 = {0.0, 1.0e-6};
  double t0 = 1.0e-6;
  // Here we use 'exact' initial points; approximations work fine too
  double f0 = 1.0; // y(t0);
  double g0 = 0.0; // dy(t0);
  std::cout << "Initial points: t0=" << t0 << ", f0=" << f0 << ", g0=" << g0
            << "\n";

  ode.solve_initial_K(t0, f0, g0);

  // Print initial points, compare to exact
  double t = t0;
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    printf("%.2f %.10f [%.10f] %.10f [%.10f]\n", t, ode.f[i], y(t), ode.g[i],
           dy(t));
    // printf("%.10f %.10f %.10f %.10f \n", t.real(), t.imag(),
    // ode.f[i].real(),
    //        ode.f[i].imag());
    t += dt;
  }

  // Solve the ODE out to very large t
  t = ode.last_t();
  std::cout << "...\n";
  for (std::size_t i = 0; i < 250; ++i) {
    t += dt;
    ode.drive();
    // if (i % 5 == 0)
    //   printf("%.10f %.10f %.10f  %.10f  \n", t.real(), t.imag(),
    //          ode.f.back().real(), ode.f.back().imag());
  }

  auto eps = (ode.f.back() - y(t)) / y(t);
  std::cout << t << " " << ode.f.back() << " " << y(t) << " " << eps << '\n';
  printf("%.2f %.10f [%.10f] %.10f [%.10f]  %.2e\n", t, ode.f.back(), y(t),
         ode.g.back(), dy(t), eps);
  // printf("%.10f %.10f %.10f  %.10f  \n", t.real(), t.imag(),
  //        ode.f.back().real(), ode.f.back().imag());
}