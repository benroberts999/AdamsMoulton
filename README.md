# Adams Moulton

A single-file header-only implementation of the Adams-Moulton method for solving Ordinary Differential Equations (ODEs) written using modern c++. Requires c++17.

Solves a general 2D system of ODEs using K-step Adams-Moutlon method,
including those with inhomogenous terms.
K is implemented from 1 to 12.
Works for real and complex value problems.

A 2D system of ODEs could be, e.g., a 2nd-order ODE, or a pair of coupled first-order ODEs.

Could reaonably simple be extened to general N-dimension problems.

## Contents

1. [Compilation and inclusion](#1-compilation-and-inclusion)
2. [Definition of problem](#2-definition-of-problem)
3. [Using the method](#3-using-the-method)
4. [Examples](#4-examples)
5. [Full documentation](#5-full-documentation)

---

## 1. Compilation and inclusion

As it's implemented as a single-file header-only library, inclusion into your project is as simple as downloading `AdamsMoulton.hpp`, and adding the include directive:

```cpp
#include "path/to/AdamsMoulton.hpp"
```

It requires c++17 to work, and has been tested with g++ and clang++.
The other files included in the repository are examples and tests, which are not required.

The code is documented using doxygen-style comments, meaning it will integrate autmatically into existing doxygen-generated documentation.

---

## 2. Definition of problem

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
F_{n+K} = F_{n+K-1} + \delta t \sum_{i=0}^K a_i \left[\frac{dF}{dt}\right]_{n+i}.
$$

Separating the $i=K$ term from the sum allows us to express $F_{n+K}$ in terms of the function values (and derivatives) at the previous $K$ points:

$$
F_{n+K} =
\left[1-\delta t \,a_K D_{n+K}\right]^{-1}
\left(F_{n+K-1} + \delta t \sum_{i=0}^{K-1} a_i \left[\frac{dF}{dt}\right]_{n+i}\right),
$$

which involves a matrix inversion (in this case, D is a 2x2 matrix).
Therefore, in a $K$-step AM method, it is required that the previous $K$ points are known.

A relatively accurate and simple way to generate the first $K$ points given just a single initial value is to use successive M-step AM methods with $M=1,2,...,K-1$. This library does that automatically.

The $a_k$ coeficients may be generated using the formula:

$$
a_k = \frac{(-1)^{K-k}}{k!(K-k)!} \int_0^1\frac{\Pi_{i=0}^K(u+i-1)}{u+K-k-1}{\rm d}u
$$

These are coded in for $K=[1,12]$, but may easily be extended.

---

## 3. Using the method

There are four steps to using the method.

  1. Define ODE system by implementing `DerivativeMatrix`
  2. Construct the solver with fixed step-size
  3. Set the initial point(s)
  4. Drive the ODE and extract the solution

Here, we give an overview to all the required functions.
It might be easier to start by following the simple example in the next setion, or looking at the examples included in this repository (in /examples/).

#### 1. Define ODE system

First, we define the ODE system.
We do this by providing an implementation to the pure-virtual struct `DerivativeMatrix`. Define your own struct which derives from `DerivativeMatrix`, and implement the a,b,c,d functions.

```cpp
  struct ODEDerivative : AdamsMoulton::DerivativeMatrix<T, Y> {
    Y a(T t) const final { return /*...*/; }
    Y b(T t) const final { return /*...*/; }
    Y c(T t) const final { return /*...*/; }
    Y d(T t) const final { return /*...*/; }
    Y Sf(T t) const final { return /*...*/; } // optional
    Y Sg(T t) const final { return /*...*/; } // optional
  };
```

* The template parameter `Y` is the return value of the D matrix; it must match the type of the ODE function F(t)
  * it is usually `double`, but may also be `float` or `complex<double>` etc.
* The template parameter `T` is the type of the _argument_ of the D matrix; it must match the type of the _argument_ of the ODE function F(t)
  * it is usually the same as `Y`, but may also be an integral/index type, e.g., when the derivative matrix is only known numerically at certain grid points/stored in an array
* The `a,b,c,d` functions _must_ be implemented. The `Sf` and `Sg` functions, which defiine the inhomogenous term, however, are optional.
* The `final` keyword is optional, but recommended, as it allows some compiler optimisations provided you never derive from `ODEDerivative`). The `const` keyword is not optional.
* Here we named the struct `ODEDerivative`, but you may name it anything

#### 2. Construct the solver

Then, we construct the solver with fixed step-size, dt.

```cpp
  ODEDerivative D{/*possible arguments*/};
  AdamsMoulton::ODESolver_2x2<K, T, Y> ode{dt, &D};
```

* `D` is an instantiation of the user-defined `DerivativeMatrix` struct
* `K` is a compile-time constant corresponding to the order of the AM method (1 to 12)
* The template parameters `Y` and `T` are the same as those from `DerivativeMatrix`
* `dt` is a constant step-size.
  * `dt` must be constant, since it must remain consistant between the K+1 and previous K points. It may, however, be positive or negative, real or complex
  * It's perfectly possible to have a non-uniform step-sizes - this introduces a Jacobian into the Derivative matrix; dt must still be constant. One of the examples shows a case like this
  * Note that `dt` has type `Y`, not `T` (though in most cases `T`=`Y`; see the examples for a case where this is not so)
* The solver, `ODESolver_2x2` takes (and stores) a pointer to an instantiation of the user-defined `DerivativeMatrix` struct. This instantiation must therefore outlive the ``ODESolver_2x2` object.

#### 3. Set the initial conditions

The first K points of the function F, and derivative dF/dt, must be known.

You may use the provided function to do this:

```cpp
  void solve_initial_K(T t0, Y f0, Y g0);
```

* t0 is the initial point
* f0 is the initial value of $f(t)$ (first component of $F(t)$)
* g0 is the initial value of $g(t)$ (second component of $F(t)$)

This automatically sets the first K values for F (and dF), given a single
initial value for F, f0=f(t0), fg=g(t0), by using successive N-step AM
methods, for N={1,2,...,K-1}.

Alternatively, you may directly access the f,g (function) and df,dg (derivative) arrays, to set these points manually:

```cpp
  for (std::size_t i = 0; i < ode.K_steps(); ++i) {
    ode.f.at(i) = /*value*/;
    ode.g.at(i) = /*value*/;
    ode.df.at(i) = /*value*/;
    ode.dg.at(i) = /*value*/;
  }
```

* `f`,`g`,`df`, and `dg` are publically-accessible arrays
* These size $K$ arrays of type $Y$ (`std::array<Y, K>`)
* `f` and `g` hold the function values at the previous $K$ points
* `df`, and `dg` hold the derivatives at the previous $K$ points
* When you drive the ODE forward, only the last $K$ points are kept. You probably want to extract these solutions and store them however you regularly would (see example)

#### 4. Drive the ODE

There are two functions to drive the ODE forwards. They do the exact same thing, which you use depends on your situation

```cpp
void drive();
void drive(T t);
```

* The first function automatically drives the ODE forards to the next `t` value
  * The next `t` value will be `t_prev+dt` (where `t_prev` is the previous `t` value) if `t` is arithetic (double, complex)
  * Or, if `t` is integer (corresponding to an array index) it will be `t+1` or `t-1`, depending on the sign of `dt`
* The second does the same, but will evaluate the derivative matrix at the specific value of `t` provided. Note: this should be used with care: it is still assumes that this value corresponds to the
  * The reson to call this is to avoid build-up or errors stemming from grids that don't align exactly.
  * For example, the 10,000th point along `t` grid may not exactly line up with `t0 * 10000*dt`, particularly for complicated non-uniform grids.

* These functions will solves the ODE at the new `t` value
* They assume system has already been solved for the K previous values {t-K*dt, ..., t-dt}.
* They also update the `f`, `g`, `df`, and `dg` arrays, which always hold the previous $K$ solutions
* The final value in each array corresponds to the most recent solution
* The order of these arrays is always the same, even if we are driving the ODE backwards

You can extract the most recent solutions using these functions:

```cpp
  // Returns most recent f value:
  Y last_f();
  // Returns most recent g value:
  Y last_g();
  // Returns most recent t value:
  T last_t();
```

You may also directly access the `f`, `g`, `df`, and `dg` arrays

---

## 4. Examples

There are several examples provided, located in the /examples/ directory.
They may be compiles using the simple provided Makefile.
Each is designed the demonstrate a capability of the library.

* **Bessel** -- [Bessel.cpp](examples/Bessel.cpp)
  * Solves the Bessel equation: a common second-order ODE
  * Demonstrates driving ODE forwards (dt>0) and backwards (dt<0)
    * This example uses GSL (GNU Scientific Library), just to compare against the expected result
    * It can be installed, e.g., on ubuntu: `apt install libgsl-dev`
* **Complex**
  * Demonstrates use of complex numbers
* **Inhomogenous**
  * Demonstrates ODE with inhomogenous term
* **Schrodinger**
  * Solves the Schrodinger for hydrogen. Demonstrates an ODE with:
    * a non-uniformly spaced grid (i.e., introduces a Jacobian)
    * and, a case where DerivativeMatrix is only known as given grid points (i.e., stored on an array)
    * In this case, while `Y=double`, `T` is an integer (array index) type (`T=std::size_t`)
* **Dirac**
  * Solves the Dirac equation for hydrogen. Demonstrates ODE with:
  * Demonstrates the case of a pair of coupled first-order ODEs (in all other examples, the system of ODEs is a single second-order ODE)

You should look into the examples in `/examples/` directory for the full example.
Here, we just show a basic outline of solving essel's differential equation

#### Example: Bessel's differential equation

* See `/examples/Bessel.cpp` for full working example.

$$
  x^2 y''(x)+ x y'(x) + (x^2-n^2)y = 0
$$

With $y(0) = 1.0$ and $y'(0)=0$, the solutions are the Bessel functions, y(x) = J_n(x)

This can be re-written as a pair of coupled first-order ODEs:

$$
  \begin{align}
  \frac{dy}{dx} &\equiv p \\
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

  // Print the first K points:
  double t = t0;
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    std::cout << t << " " << ode.f.at(i) << '\n';
    t += ode.dt();
  }

  // Drive forwards another 100 steps
  for (int i = 0; i < 100; ++i) {
    ode.drive();
    std::cout << ode.last_t() << " " << ode.last_f() << '\n';
  }
```

---

## 5. Full Documentation

```cpp
namespace AdamsMoulton {}
```

* Contains all the classes and functions which use general N-step Adams Moulton method to solve systems of 2x2 ODEs, up to N=12.
* Everything is in this namespace

### Misc. functions/constants

```cpp
template <typename T> constexpr bool is_complex_v = is_complex<T>::value;
```

* User-defined type-trait: Checks whether T is a std::complex type
* For example:

  ```cpp
    static_assert(!is_complex_v<double>);
    static_assert(!is_complex_v<float>);
    static_assert(is_complex_v<std::complex<double>>);
    static_assert(is_complex_v<std::complex<float>>);
  ```

```cpp
template <typename T, typename U, std::size_t N, std::size_t M>
constexpr T inner_product(const std::array<T, N> &a, const std::array<U, M> &b);
```

* Inner product of two std::arrays.

$$ \mathrm{inner\_product}(a, b) = \sum_{i=0}^{N-1} a_i * b_i $$

* Where `N = min(a.size(), b.size())`.
* The types of the arrays may be different (T and U).
* However, U must be convertable to T; the return-type is T (same as first array).

---

### Adams-Moulton coefficients

```cpp
static constexpr std::size_t K_max;
```

* Stores maximum K (order of AM method) for which we have coefficients implemented (currently 12).

```cpp
template <std::size_t K>
struct AM_Coefs{};
```

* Holds the K+1 Adams-Moulton ak coefficients for the K-step AM method. Final one, aK, is stored separately.
* `K` is a compile-time constant corresponding to order of method
  * Requires `K>0` and `K<AdamsMoulton::K_max`
* The Adams coefficients, a_k, are defined such that:

 $$ F_{n+K} = F_{n+K-1} + dt * \sum_{k=0}^K a_k y_{n+k} $$

 $$ y = dF(t)/dt $$

```cpp
  static constexpr std::array<double, K> AM_Coefs::ak;
  static constexpr double AM_Coefs::aK;
```

* `ak` stores the first K coefficients: ak for k={0,1,...,K-1}
* `aK` stores the final aK coefficients: ak for k=K
* They are stored as `double` regardless of other template parameters.
* Note: the 'order' of the coefs is reversed compared to some sources.

---

### The derivative matrix (defining the ODE)

```cpp
template <typename T = double, typename Y = double> struct DerivativeMatrix{};
```

* Pure-virtual struct, holds the derivative matrix for 2x2 system of ODEs.
* Derive from this, and implement a(t),b(t),c(t),d(t) to define the 2x2 ODE.
* The template parameter `Y` is the return value of the D matrix; it must match the type of the ODE function F(t)
  * it is usually `double`, but may also be `float` or `complex<double>` etc.
* The template parameter `T` is the type of the _argument_ of the D matrix; it must match the type of the _argument_ of the ODE function F(t)
  * it is usually the same as `Y`, but may also be an integral/index type, e.g., when the derivative matrix is only known numerically at certain grid points/stored in an array

```cpp
  virtual DerivativeMatrix::Y a(T t) const;
  virtual DerivativeMatrix::Y b(T t) const;
  virtual DerivativeMatrix::Y c(T t) const;
  virtual DerivativeMatrix::Y d(T t) const;
```

* These functions define the derivative matrix D (see definition above)
* The `a,b,c,d` functions are pure virtual and _must_ be implemented.

```cpp
  virtual Y DerivativeMatrix::Sf(T) const;
  virtual Y DerivativeMatrix::Sg(T) constk
```

* The `Sf` and `Sg` functions, define the inhomogenous term.
* These may optionally be over-written in the derived struct to add an inhomogenous term.
* These are optional (be default, they return 0).

---

### The ODE solver

```cpp
template <std::size_t K, typename T = double, typename Y = double>
class ODESolver_2x2{};
```

* Solves a 2x2 system of ODEs using a K-step Adams-Moutlon method
* Form of the ODE defined above
* Template parameter, `T`, is the type of the argument of the Derivative Matrix (i.e., type of `t`).  
  * This is often `double` or `complex<double>`, but may also be an index type (e.g., std::size_t) if the derivative matrix is only known numerically at certain grid points/stored in an array.
* Template parameter, `Y`, is the type of the function value F(t), and the
type of dt, and the return value of the Derivative Matrix.
  * This is often `double`, but may also be another floating-point type, or std::complex.

#### Constructor

```cpp
  ODESolver_2x2(Y dt, const DerivativeMatrix<T, Y> *D);
```

* Constructor
* Takes in step-size, `dt`
  * Type `Y` (not `T`)
  * Is a constant value; may not change inside the solver
  * The step-size, dt, may be positive (to drive forwards) or negative (to drive backwards); it may also be complex.
* Takes in a pointer to the DerivativeMatrix
  * This pointer to the DerivativeMatrix is stored.
  * This may not be null,and must outlive the `ODESolver_2x2`.

#### Public data members

```cpp
  std::array<Y, K> ODESolver_2x2::f;
  std::array<Y, K> ODESolver_2x2::g;
  std::array<Y, K> ODESolver_2x2::df;
  std::array<Y, K> ODESolver_2x2::dg;
```

* Arrays of size K, type Y.
* f and g hold the previous K functions values
* g and f hold the previous K derivative values
* Each value is separated by dt

#### Public member functions

```cpp
  constexpr std::size_t ODESolver_2x2::K_steps() const;
```

* Returns the AM order (number of steps), K

```cpp
Y ODESolver_2x2::last_f();
Y ODESolver_2x2::last_g();
```

* Returns most recent values for f and g.
* You may also access the f/g arrays directly

```cpp
T ODESolver_2x2::last_t();
```

* Returns most recent t value
* nb: last_f() := f(last_t())

```cpp
Y ODESolver_2x2::dt() { return m_dt; }
```

* Returns the step size

```cpp
Y ODESolver_2x2::dfdt(Y ft, Y gt, T t) const;
Y ODESolver_2x2::dgdt(Y ft, Y gt, T t) const;
```

* Returns derivative, df/dt(t) and df/dt(t), given f(t),g(t),t
* Invokes the DerivativeMatrix

```cpp
void ODESolver_2x2::drive();
void ODESolver_2x2::drive(T t);
```

* Drives the DE system to next value, F(t), assuming system has already been solved for the K previous values {t-K*dt, ..., t-dt}.
* Note: t must align properly with previous values: i.e., t = last_t + dt \n
  * We may re-send t to the function to avoid large build-up of small errors. For example, the 10,000th point along t grid may not exactly line up with `t0 * 10000*dt`, particularly for complicated non-linear grids.
  * There's an overload that avoids this, but care should be taken.
* The type of t (`T`) must match type required to compute DerivativeMatrix.
  