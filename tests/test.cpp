#define CATCH_CONFIG_MAIN
#include "AdamsMoulton.hpp"
#include "tests/catch2/catch.hpp"
#include <cmath>

//------------------------------------------------------------------------------
TEST_CASE("Simple: integrate forwards") {

  // y''(x) = -y(x)
  // y(0)=0, y'(0)=1
  // => y(x) = sin(x)
  struct DM : AdamsMoulton::DerivativeMatrix<double> {
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double) const final { return -1.0; }
    double d(double) const final { return 0.0; }
  };
  DM D;

  double dt = 0.001;
  AdamsMoulton::ODESolver_2x2<5> ode{dt, &D};

  REQUIRE(ode.dt() == dt);

  double t0 = 0.0;
  double f0 = 0.0;
  double g0 = 1.0;
  ode.solve_initial_K(t0, f0, g0);

  // test the initial points:
  double t = 0.0;
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    REQUIRE(ode.f[i] == Approx(std::sin(t)).epsilon(1.0e-5));
    REQUIRE(ode.g[i] == Approx(std::cos(t)).epsilon(1.0e-5));
    t += ode.dt();
  }
  REQUIRE(ode.last_f() == ode.f.back());
  REQUIRE(ode.last_g() == ode.g.back());

  int num_steps = 1000;
  for (int i = 0; i < num_steps; ++i) {
    ode.drive();
  }
  REQUIRE(ode.last_f() == Approx(std::sin(ode.last_t())).epsilon(1.0e-5));
  REQUIRE(ode.last_g() == Approx(std::cos(ode.last_t())).epsilon(1.0e-5));

  REQUIRE(ode.last_t() ==
          Approx(t0 + (num_steps + (int)ode.K_steps() - 1) * ode.dt()));
}

//------------------------------------------------------------------------------
TEST_CASE("Simple: integrate backwards") {

  // y''(x) = -y(x)
  // y(0)=0, y'(0)=1
  // => y(x) = sin(x)
  struct DM : AdamsMoulton::DerivativeMatrix<double> {
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double) const final { return -1.0; }
    double d(double) const final { return 0.0; }
  };
  DM D;

  double dt = -0.001;
  AdamsMoulton::ODESolver_2x2<5> ode{dt, &D};

  double t0 = 0.0;
  double f0 = 0.0;
  double g0 = 1.0;
  ode.solve_initial_K(t0, f0, g0);

  // test the initial points:
  double t = 0.0;
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    REQUIRE(ode.f[i] == Approx(std::sin(t)).epsilon(1.0e-5));
    REQUIRE(ode.g[i] == Approx(std::cos(t)).epsilon(1.0e-5));
    t += ode.dt();
  }

  int num_steps = 1000;
  for (int i = 0; i < num_steps; ++i) {
    ode.drive();
  }
  REQUIRE(ode.last_f() == Approx(std::sin(ode.last_t())).epsilon(1.0e-5));
  REQUIRE(ode.last_g() == Approx(std::cos(ode.last_t())).epsilon(1.0e-5));
  REQUIRE(ode.last_t() < t0);

  REQUIRE(ode.last_t() ==
          Approx(t0 + (num_steps + (int)ode.K_steps() - 1) * ode.dt()));
}

//------------------------------------------------------------------------------
TEST_CASE("Simple: use drive(t)") {

  // y''(x) = -y(x)
  // y(0)=0, y'(0)=1
  // => y(x) = sin(x)
  struct DM : AdamsMoulton::DerivativeMatrix<double> {
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double) const final { return -1.0; }
    double d(double) const final { return 0.0; }
  };
  DM D;

  double dt = 0.001;
  AdamsMoulton::ODESolver_2x2<5> ode{dt, &D};

  REQUIRE(ode.dt() == dt);

  double t0 = 0.0;
  double f0 = 0.0;
  double g0 = 1.0;
  ode.solve_initial_K(t0, f0, g0);

  // test the initial points:
  double t = 0.0;
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    REQUIRE(ode.f[i] == Approx(std::sin(t)).epsilon(1.0e-5));
    REQUIRE(ode.g[i] == Approx(std::cos(t)).epsilon(1.0e-5));
    t += ode.dt();
  }
  REQUIRE(ode.last_f() == ode.f.back());
  REQUIRE(ode.last_g() == ode.g.back());

  int num_steps = 1000;
  for (int i = 0; i < num_steps; ++i) {
    ode.drive(t);
    t += ode.dt();
  }
  REQUIRE(ode.last_f() == Approx(std::sin(ode.last_t())).epsilon(1.0e-5));
  REQUIRE(ode.last_g() == Approx(std::cos(ode.last_t())).epsilon(1.0e-5));

  REQUIRE(ode.last_t() ==
          Approx(t0 + (num_steps + (int)ode.K_steps() - 1) * ode.dt()));
  REQUIRE(ode.last_t() == t - ode.dt());
}

//------------------------------------------------------------------------------
TEST_CASE("Simple: array index argument") {

  // y''(x) = -y(x)
  // y(0)=0, y'(0)=1
  // => y(x) = sin(x)
  struct DM : AdamsMoulton::DerivativeMatrix<int, double> {
    double a(int) const final { return 0.0; }
    double b(int) const final { return 1.0; }
    double c(int) const final { return -1.0; }
    double d(int) const final { return 0.0; }
  };
  DM D;

  double dt = 0.001;
  AdamsMoulton::ODESolver_2x2<5, int, double> ode{dt, &D};

  REQUIRE(ode.dt() == dt);

  int t0 = 0;
  double f0 = 0.0;
  double g0 = 1.0;
  ode.solve_initial_K(t0, f0, g0);

  // test the initial points:
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    double t = double(i) * ode.dt();
    REQUIRE(ode.f[i] == Approx(std::sin(t)).epsilon(1.0e-5));
    REQUIRE(ode.g[i] == Approx(std::cos(t)).epsilon(1.0e-5));
  }
  REQUIRE(ode.last_f() == ode.f.back());
  REQUIRE(ode.last_g() == ode.g.back());

  int num_steps = 1000;
  for (int i = 0; i < num_steps; ++i) {
    ode.drive();
  }
  REQUIRE(ode.last_f() ==
          Approx(std::sin(ode.last_t() * ode.dt())).epsilon(1.0e-5));
  REQUIRE(ode.last_g() ==
          Approx(std::cos(ode.last_t() * ode.dt())).epsilon(1.0e-5));

  REQUIRE(ode.last_t() == Approx((num_steps + (int)ode.K_steps() - 1)));
}

//------------------------------------------------------------------------------
TEST_CASE("Simple: K=1,12") {
  // y''(x) = -y(x)
  // y(0)=0, y'(0)=1
  // => y(x) = sin(x)
  struct DM : AdamsMoulton::DerivativeMatrix<double> {
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double) const final { return -1.0; }
    double d(double) const final { return 0.0; }
  };
  DM D;

  double dt = 0.001;
  AdamsMoulton::ODESolver_2x2<1> ode_1{dt, &D};
  AdamsMoulton::ODESolver_2x2<12> ode_12{dt, &D};

  double t0 = 0.0;
  double f0 = 0.0;
  double g0 = 1.0;
  ode_1.solve_initial_K(t0, f0, g0);
  ode_12.solve_initial_K(t0, f0, g0);

  // test the initial points:
  double t = 0.0;
  for (std::size_t i = 0; i < ode_1.f.size(); ++i) {
    REQUIRE(ode_1.f[i] == Approx(std::sin(t)).epsilon(1.0e-4));
    REQUIRE(ode_1.g[i] == Approx(std::cos(t)).epsilon(1.0e-4));
    t += ode_1.dt();
  }
  t = 0.0;
  for (std::size_t i = 0; i < ode_12.f.size(); ++i) {
    REQUIRE(ode_12.f[i] == Approx(std::sin(t)).epsilon(1.0e-5));
    REQUIRE(ode_12.g[i] == Approx(std::cos(t)).epsilon(1.0e-5));
    t += ode_12.dt();
  }

  int num_steps = 1000;
  for (int i = 0; i < num_steps; ++i) {
    ode_1.drive();
    ode_12.drive();
  }
  REQUIRE(ode_1.last_f() == Approx(std::sin(ode_1.last_t())).epsilon(1.0e-4));
  REQUIRE(ode_12.last_f() == Approx(std::sin(ode_12.last_t())).epsilon(1.0e-5));
}

//------------------------------------------------------------------------------
TEST_CASE("Complex") {

  /*
  y''(z) = - i*y(z)
  The exact solution, with y(0)=0, y'(0)=1, is:
  y(z) = - (-1)^(3/4) * sin[(-1)^(1/4)*z]
  */

  using namespace std::complex_literals;

  auto y_expected = [](auto x) {
    const auto r2 = std::sqrt(2.0);
    const auto pow1 = -1.0 / r2 + 1.0i / r2; // (-1)^(3/4)
    const auto pow2 = 1.0 / r2 + 1.0i / r2;  // (-1)^(1/4)
    return -pow1 * std::sin(pow2 * x);
  };

  // Define the DerivativeMatrix for the Bessel equation
  struct ComplexDerivative
      : AdamsMoulton::DerivativeMatrix<std::complex<double>,
                                       std::complex<double>> {
    std::complex<double> a(std::complex<double>) const final { return 0.0; }
    std::complex<double> b(std::complex<double>) const final { return 1.0; }
    std::complex<double> c(std::complex<double>) const final { return -1.0i; }
    std::complex<double> d(std::complex<double>) const final { return 0.0; }
  };

  ComplexDerivative D{};

  constexpr std::size_t N_step = 5; // use 5-step method

  using ComplexAdams = AdamsMoulton::ODESolver_2x2<N_step, std::complex<double>,
                                                   std::complex<double>>;

  std::complex<double> z0 = 0.0;
  std::complex<double> y0 = 0.0;
  std::complex<double> dy0 = 1.0;

  // First, solve along the positive real axis, starting at 0:
  std::complex<double> dz = {1.0e-3, 0.0};
  auto ode = ComplexAdams{dz, &D};
  ode.solve_initial_K(z0, y0, dy0);

  // Then, solve along the positive imaginary axis, starting at 0:
  std::complex<double> idz = {0.0, 1.0e-3};
  auto ode_iz = ComplexAdams{idz, &D};
  ode_iz.solve_initial_K(z0, y0, dy0);

  // Integrate outwards to |z|=10
  double abs_z_max = 10.0;
  while (std::abs(ode.last_t()) < abs_z_max) {
    ode.drive();
    ode_iz.drive();
  }

  REQUIRE(std::abs(ode.last_t()) ==
          Approx(abs_z_max).margin(1.01 * std::abs(dz)));

  REQUIRE(ode.last_f().real() ==
          Approx(y_expected(ode.last_t()).real()).epsilon(1.0e-5));
  REQUIRE(ode.last_f().imag() ==
          Approx(y_expected(ode.last_t()).imag()).epsilon(1.0e-5));

  REQUIRE(ode_iz.last_f().real() ==
          Approx(y_expected(ode_iz.last_t()).real()).epsilon(1.0e-5));
  REQUIRE(ode_iz.last_f().imag() ==
          Approx(y_expected(ode_iz.last_t()).imag()).epsilon(1.0e-5));
}