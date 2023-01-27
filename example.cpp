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

  struct DD2 : AdamsMoulton::DerivativeMatrix<double> {
    double a(double) const { return 0.0; }
    double b(double) const { return 1.0; }
    double c(double) const { return -1.0; }
    double d(double t) const { return -1.0 / t; }
  };

  DD<double> D1{};
  DD2 D2{};

  double t = 10.0;
  double dt = 0.01;

  AdamsMoulton::ODESolver_2x2<12> a{dt, &D2};

  double f0 = y(t);
  double g0 = dy(t);
  a.initial(t, f0, g0);

  for (std::size_t i = 0; i < a.f.size(); ++i) {
    // a.f[i] = y(t);
    // a.g[i] = dy(t);
    // a.df[i] = a.dfdt(a.f[i], a.g[i], t);
    // a.dg[i] = a.dgdt(a.f[i], a.g[i], t);
    // std::cout << t << " " << a.f[i] << " " << y(t) << '\n';
    printf("%.2f %.10f [%.10f] %.10f [%.10f]\n", t, a.f[i], y(t), a.g[i],
           dy(t));
    t += dt;
  }

  std::cout << "--\n";
  for (std::size_t i = 0; i < 1000000; ++i) {
    a.drive(t);

    // auto eps = (a.f.back() - y(t)) / y(t);
    // std::cout << t << " " << a.f.back() << " " << y(t) << " " << eps << '\n';

    t += dt;
  }

  t -= dt;
  auto eps = (a.f.back() - y(t)) / y(t);
  // std::cout << t << " " << a.f.back() << " " << y(t) << " " << eps << '\n';
  printf("%.2f %.10f [%.10f] %.10f [%.10f]  %.2e\n", t, a.f.back(), y(t),
         a.g.back(), dy(t), eps);
}