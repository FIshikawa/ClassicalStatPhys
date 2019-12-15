#ifndef INTEGRATOR_EULER_HPP
#define INTEGRATOR_EULER_HPP

#include <string>
#include <vector>

namespace integrator {

class Euler {
public:
  static std::string name() { return "Euler method"; }
  Euler(unsigned int dim) : dim_(dim), k_(dim) {}
  template<class F>
  void step(double t, double h, std::vector<double>& y, F const& f) const {
    f(t, y, k_);
    for (int i = 0; i < dim_; ++i) y[i] += h * k_[i];
  }
private:
  unsigned int dim_;
  mutable std::vector<double> k_;
};

} // namespace integrator

#endif // INTEGRATOR_EULER_HPP
