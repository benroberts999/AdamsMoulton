// Simple c++ class for doing approximate comparisons (for testing)
#pragma once
#include <cmath>

namespace qip {

template <typename T> class Approx {
private:
  struct Approx_EPS {
    T value;
    T epsilon;
    friend bool operator==(T t, Approx_EPS e) {
      return std::abs((t - e.value) / e.value) <= e.epsilon;
    }
  };

  struct Approx_DEL {
    T value;
    T delta;
    friend bool operator==(T t, Approx_DEL d) {
      return std::abs(t - d.value) <= d.delta;
    }
  };

public:
  T value;
  Approx(T t) : value(t) {}
  Approx_EPS eps(T epsilon) const { return Approx_EPS{value, epsilon}; }
  Approx_DEL del(T delta) const { return Approx_DEL{value, delta}; }

  friend bool operator==(T t, Approx<T> a) {
    return std::abs(t - a.value) <= 1.0e-14;
  }
};

} // namespace qip
