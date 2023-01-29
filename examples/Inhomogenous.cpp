#include "AdamsMoulton.hpp"
#include <cmath>
#include <complex>
#include <iostream>

// Note: This example uses GSL (GNU Scientific Libraries) for evaluating `exact'
// bessel functions; only required for testing the result.

const std::string description{R"(
Adams Moulton example: Inhomogenous ODE.

y''(x) = -y(x) + sin(x)

Can be written as:
    dF/dx = D * F + S(x)

where:
    F(x) = ( y(x)  )
           ( y'(x) ),

    D(t) = ( 0    1 )
           ( -1   0 ),
and 
    S(x) = ( sin(x) )
           ( 0      ).

The exact solution, with y(0)=0 and y'(0)=1, is: 
y(x) = 0.5 * [ 3*sin(x) - x*cos(x) ]
)"};

//==============================================================================
int main() {

  std::cout << description << "\n";

  // Define the DerivativeMatrix for the Bessel equation
  struct InhomogDerivative : AdamsMoulton::DerivativeMatrix<double> {
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double) const final { return -1.0; }
    double d(double) const final { return 0.0; }
    double Sf(double) const final { return 0.0; }
    double Sg(double x) const final { return std::sin(x); }
  };

  InhomogDerivative D{};

  // Exact solution: just to compare result
  auto y = [=](double x) {
    return 0.5 * (3.0 * std::sin(x) - x * std::cos(x));
  };
  auto dy = [=](double x) {
    return 0.5 * (x * std::sin(x) + 2.0 * std::cos(x));
  };

  // Set the stp size
  double dt = 0.0001;
  // Construct the solver:
  AdamsMoulton::ODESolver_2x2<5> ode{dt, &D};

  // Set the initial conditions:
  double t0 = 0.0;
  double f0 = 0.0;
  double g0 = 1.0;
  std::cout << "Initial points: t0=" << t0 << ", f0=" << f0 << ", g0=" << g0
            << "\n";
  std::cout << "Using a K=" << ode.K_steps() << " AM method\n";

  ode.solve_initial_K(t0, f0, g0);

  // Print initial points, compare to exact

  std::cout << "Compare initial K points to expected (exact) solution:\n";
  std::cout << "  t       y(t)         [Exact Soln   ]  dy/dt(t)     [Exact "
               "Soln   ]\n";
  double t = t0;
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    printf("%8.4f %13.10f [%13.10f] %13.10f [%13.10f]\n", t, ode.f[i], y(t),
           ode.g[i], dy(t));
    t += ode.dt();
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
}