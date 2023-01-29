#include "AdamsMoulton.hpp"
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

const std::string description{R"(
Adams Moulton example: complex values.

The equation:

y''(z) = - i*y(z)

Can be written as:
    dF/dz = D * F

where:
    F(z) = ( y(z)  )
           ( dy/dz )
and:
    D(z) = ( 0   1 )
           (-i   0 ).

The exact solution, with y(0)=0, y'(0)=1, is:
y(z) = - (-1)^(3/4) * sin[(-1)^(1/4)*z]
)"};

//==============================================================================
int main() {

  std::cout << description << "\n";

  using namespace std::complex_literals;

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

  constexpr std::size_t N_step = 12; // use 12-step method

  using ComplexAdams = AdamsMoulton::ODESolver_2x2<N_step, std::complex<double>,
                                                   std::complex<double>>;

  // We will solve the ODE twice, once along the real axis, then along the
  // imaginary axis.
  // Both times, start from (0.0 + 0.0i), using y(0)=0, y'(0)=1 + 0i
  std::complex<double> z0 = 0.0;
  std::complex<double> y0 = 0.0;
  std::complex<double> dy0 = 1.0;

  // vectors to store the solutions. Along the real and imaginary axis
  std::vector<std::pair<double, std::complex<double>>> y_out_realz;
  std::vector<std::pair<double, std::complex<double>>> y_out_imagz;

  // First, solve along the positive real axis, starting at 0:
  std::complex<double> dz = {1.0e-3, 0.0};
  auto ode = ComplexAdams{dz, &D};
  ode.solve_initial_K(z0, y0, dy0);

  // Then, solve along the positive imaginary axis, starting at 0:
  std::complex<double> idz = {0.0, 1.0e-3};
  auto ode_iz = ComplexAdams{idz, &D};
  ode_iz.solve_initial_K(z0, y0, dy0);

  // store the initial K points in the output array:
  auto rz = z0;
  auto iz = z0;
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    y_out_realz.emplace_back(rz.real(), ode.f.at(i));
    y_out_imagz.emplace_back(iz.imag(), ode_iz.f.at(i));
    rz += ode.dt();
    iz += ode_iz.dt();
  }

  // Integrate outwards to |z|=10
  double abs_z_max = 10.0;
  std::cout << "Solving out to |z|=" << abs_z_max
            << ", along thr real and imaginary axis\n";
  while (std::abs(ode.last_t()) < abs_z_max) {
    ode.drive();
    ode_iz.drive();
    y_out_realz.emplace_back(ode.last_t().real(), ode.last_f());
    y_out_imagz.emplace_back(ode_iz.last_t().imag(), ode_iz.last_f());
  }

  // Output result to files for plotting:

  std::string re_fname{"y_real_axis.txt"};
  std::string im_fname{"y_imag_axis.txt"};
  std::ofstream re_out(re_fname);
  std::ofstream im_out(im_fname);

  auto y_expected = [](auto x) {
    const auto r2 = std::sqrt(2.0);
    const auto pow1 = -1.0 / r2 + 1.0i / r2; // (-1)^(3/4)
    const auto pow2 = 1.0 / r2 + 1.0i / r2;  // (-1)^(1/4)
    return -pow1 * std::sin(pow2 * x);
  };

  re_out << "Re(z)     Re(y)     Im(y)     Re(y_exact)     Im(y_exact)\n";
  for (auto &[z, y] : y_out_realz) {
    auto y_exact = y_expected(z);
    re_out << z << " " << y.real() << " " << y.imag() << " " << y_exact.real()
           << " " << y_exact.imag() << "\n";
  }

  im_out << "Im(z)     Re(y)     Im(y)     Re(y_exact)     Im(y_exact)\n";
  for (auto &[z, y] : y_out_imagz) {
    auto y_exact = y_expected(1.0i * z);
    im_out << z << " " << y.real() << " " << y.imag() << " " << y_exact.real()
           << " " << y_exact.imag() << "\n";
  }

  std::cout << "\nSolutions written to files:\n"
            << re_fname << " - for real axis, and \n"
            << im_fname << " - for imag. axis.\n";

  std::cout << "These may be plotted, for example, using gnuplot with:\n";

  const std::string gnuplot{R"(
   set title 'y(z), with Im(z)=0'
   set xlabel 'Re(z)'
   plot for [i=2:5] 'y_real_axis.txt' u 1:i w l t columnheader(i)

   set title 'y(z), with Re(z)=0'
   set xlabel 'Im(z)'
   plot for [i=2:5] 'y_imag_axis.txt' u 1:i w l t columnheader(i)
)"};

  std::cout << gnuplot << "\n";
}
