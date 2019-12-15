#ifndef INTEGRATOR_RUNGE_KUTTA_2ND_HPP
#define INTEGRATOR_RUNGE_KUTTA_2ND_HPP

#include <string>
#include <vector>

namespace integrator {

class RungeKutta2nd{
public:
  static std::string name() { return "2nd-order Runge-Kutta method"; }
  RungeKutta2nd(unsigned int dim) : dim_(dim), k1_(dim), k2_(dim), yt_(dim) {}
  template<class F>
  void step(double t, double h, std::vector<double>& y, F const& f) const {
    double h2 = h / 2;
    f(t, y, k1_);
    for (int i = 0; i < dim_; ++i) yt_[i] = y[i] + h2 * k1_[i];
    f(t + h2, yt_, k2_);
    for (int i = 0; i < dim_; ++i) y[i] += h * k2_[i];
  }
private:
  unsigned int dim_;
  mutable std::vector<double> k1_;
  mutable std::vector<double> k2_;
  mutable std::vector<double> yt_;
};
  
} // namespace integrator

#endif // INTEGRATOR_RUNGE_KUTTA_2ND_HPP
