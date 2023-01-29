# Adams Moulton

Solves a general 2D system of ODEs using K-step Adams-Moutlon method,
including those with inhomogenous terms.
K is implemented from 1 to 12.
Works for complex and real value problems.

Could reaonably simple be extened to general N-dimension problems.

## Definition of problem

The system of ODEs is defined such that:

$$ \frac{dF(t)}{dt} = D(t) * F(t) + S(t) $$

Where F is a 2D set of functions:

$$
  F(t) = \begin{pmatrix}
    f(t)\\
    g(t)
  \end{pmatrix},
$$

D is the 2x2 "derivative matrix":

$$
  D(t) = \begin{pmatrix}
    a(t) & b(t)\\
    c(t) & d(t)
  \end{pmatrix},
$$

and S(t) is the (optional) 2D inhomogenous term:

$$
  S(t) = \begin{pmatrix}
    s_f(t)\\
    s_g(t)
  \end{pmatrix}.
$$

In the Adams-Moutlon method, the ODE is written
$$
F_{n+K} = F_{n+K-1} + \delta t \sum_{i=0}^K a_i \left(\frac{dF}{dt}\right)_{n+i}.
$$
Separating the $i=K$ term from the sum allows us to express $F_{n+K}$ in terms of the function values (and derivatives) at the previous $K$ points:
$$
F_{n+K} = 
\left[1-\delta t a_K D_{n+K}\right]^{-1}
\left(F_{n+K-1} + \delta t \sum_{i=0}^{K-1} a_i \left(\frac{dF}{dt}\right)_{n+i}\right),
$$
which involves a matrix inversion (in this case, D is a 2x2 matrix).
Therefore, in a $K$-step AM method, it is required that the previous $K$ points are known.

## Using the method

1. Define ODE system by implementing `DerivativeMatrix`
2. Construct the solver with fixed step-size
3. Set the initial point(s)
4. Drive the ODE

See struct `DerivativeMatrix` - which is a pure virtual struct that must be
implmented by the user in order to define the ODE.

The step-size, dt, must remain constant (since it must remain consistant
between the K+1 and previous K points). It may be positive or negative,
however (or complex). It's perfectly possible to have a non-uniform grid - this
just introduces a Jacobian into the Derivative matrix; dt must still be
constant.

The template parameter, T, is the type of the argument of the Derivative
Matrix (i.e., type of `t`). This is often `double` ro `complex<double>`, but may
also be an index type (e.g., std::size_t) if the derivative matrix is only known
numerically at certain grid points/stored in an array.

The template parameter, Y, is the type of the function value F(t), and the
type of dt, and the return value of the Derivative Matrix. This is often
`double`, but may also be another floating-point type, or std::complex.

The first K points of the function F, and derivative dF/dt, must be known.
You may directly access the f,g (function) and df,dg (derivative) arrays, to
set these points.

Alternatively, you may use the provided function
  ```cpp
  void solve_initial_K(T t0, Y f0, Y g0);
  ```
which automatically sets the first K values for F (and dF), given a single
initial value for F, f0=f(t0), fg=g(t0), by using successive N-step AM
methods, for N={1,2,...,K-1}.

For now, just a 2x2 system. In theory, simple to generalise to an N*N system,
though requires a matrix inversion.

\par

**Example:** Bessel's differential equation

$$
  x^2 \frac{d^2y}{dx^2} + x \frac{dy}{dx} + (x^2-n^2)y = 0
$$

With y(0) = 1.0, the solutions are the Bessel functions, y(x) = J_n(x)

This can be re-written as a pair of coupled first-order ODEs:

$$
  \begin{align}
  \frac{dy}{dx} &= p \\
  \frac{dp}{dx} &= \left[\left(\frac{n}{x}\right)^2 - 1\right]y - \frac{1}{x}
p \end{align} $$

Putting this into the notation/form required for the solver we have:

$$
  F(x) = \begin{pmatrix}
    y(x)\\
    p(x)
  \end{pmatrix}
$$

with the "derivative matrix":

$$
  D(x) = \begin{pmatrix}
    0 & 1 \\
    \left(\frac{n}{x}\right)^2 - 1 & \frac{-1}{x}
  \end{pmatrix}
$$

i.e.,

for n=1:

```cpp
  struct BesselDerivative : AdamsMoulton::DerivativeMatrix<double, double> {
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double t) const final { return 1.0/t/t - 1.0; }
    double d(double t) const final { return -1.0 / t; }
  };
```

Or, more generally (for example):

```cpp
  struct BesselDerivative : AdamsMoulton::DerivativeMatrix<double, double> {
    int n;
    BesselDerivative(int tn) : n(tn) {}
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double t) const final { return std::pow(n/t,2) - 1.0; }
    double d(double t) const final { return -1.0 / t; }
  };
```

Minimal example: -- see full examples included elsewhere

```cpp
  // Construct the Derivative matrix (BesselDerivative defined above) with n=0
  int n = 0;
  BesselDerivative D{n};

  // Set the step size:
  double dt = 0.01;

  // Construct the Solver, using K=6-step method:
  AdamsMoulton::ODESolver_2x2<6> ode{dt, &D};

  // Since 1/t appears in D, we cannot start at zero. Instead, begin at small t
  double t0 = 1.0e-6;

  // Set initial points:
  // Note: these are *approximate*, since they are technically f(0.0)
  double f0 = 1.0;
  double g0 = 0.0;

  // Use automatic solver for first K points:
  ode.solve_initial_K(t0, f0, g0);

  // Drive forwards another 100 steps
  for (int i = 0; i < 100; ++i) {
    ode.drive();
    std::cout << ode.last_t() << " " << ode.last_f() << '\n';
  }
```