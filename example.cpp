#include "AdamsMoulton.hpp"
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <iostream>

double a(double) { return 0.0; }
double b(double) { return 1.0; }
double c(double) { return -1.0; }
double d(double t) { return -1.0 / t; }

template <typename T> struct DD : AdamsMoulton::DerivativeMatrix<T> {
  double a(T t) const { return ::a(t); }
  double b(T t) const { return ::b(t); }
  double c(T t) const { return ::c(t); }
  double d(T t) const { return ::d(t); }
};

//==============================================================================
int main() {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Adams Moulton test/example\n";

  int n = 0;

  auto y = [=](double x) { return gsl_sf_bessel_Jn(n, x); };

  auto dy = [=](double x) {
    return 0.5 * (gsl_sf_bessel_Jn(n - 1, x) - gsl_sf_bessel_Jn(n + 1, x));
  };

  struct BesselDerivative : AdamsMoulton::DerivativeMatrix<double> {
    int n;
    BesselDerivative(int tn) : n(tn) {}
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double t) const final { return std::pow(n / t, 2) - 1.0; }
    double d(double t) const final { return -1.0 / t; }
  };

  BesselDerivative D2{n};

  double dt = 0.01;

  AdamsMoulton::ODESolver_2x2<12> ode{dt, &D2};

  double t0 = 1.0e-4;
  // Here we use 'exact' initial points; approximations work fine too
  double f0 = y(t0);
  double g0 = dy(t0);
  std::cout << "Initial points: t0=" << t0 << ", f0=" << f0 << ", g0=" << g0
            << "\n";

  ode.solve_initial_K(t0, f0, g0);

  // Print initial points, compare to exact
  double t = t0;
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    printf("%.2f %.10f [%.10f] %.10f [%.10f]\n", t, ode.f[i], y(t), ode.g[i],
           dy(t));
    t += dt;
  }

  // Solve the ODE out to very large t
  t = ode.last_t();
  std::cout << "...\n";
  for (std::size_t i = 0; i < 1000000; ++i) {
    t += dt;
    ode.drive();
  }

  auto eps = (ode.f.back() - y(t)) / y(t);
  // std::cout << t << " " << ode.f.back() << " " << y(t) << " " << eps << '\n';
  printf("%.2f %.10f [%.10f] %.10f [%.10f]  %.2e\n", t, ode.f.back(), y(t),
         ode.g.back(), dy(t), eps);
}