#ifndef INTEGRATOR_RUNGE_KUTTA_3RD_HPP
#define INTEGRATOR_RUNGE_KUTTA_3RD_HPP

#include <string>
#include <vector>

namespace integrator {

class RungeKutta3rd{
public:
  static std::string name() { return "3rd-order Runge-Kutta method"; }
  RungeKutta3rd(unsigned int dim) : dim_(dim), k1_(dim), k2_(dim), k3_(dim), yt_(dim) {}
  template<class F>
  void step(double t, double h, std::vector<double>& y, F const& f) const {
    double h23 = 2 * h / 3;
    f(t, y, k1_);
    for (int i = 0; i < dim_; ++i) yt_[i] = y[i] + h23 * k1_[i];
    f(t + h23, yt_, k2_);
    for (int i = 0; i < dim_; ++i) yt_[i] = y[i] + h23 * k2_[i];
    f(t + h23, yt_, k3_);
    for (int i = 0; i < dim_; ++i)
      y[i] = y[i] + (h / 8) * (2 * k1_[i] + 3 * k2_[i] + 3 * k3_[i]);
  }
private:
  unsigned int dim_;
  mutable std::vector<double> k1_;
  mutable std::vector<double> k2_;
  mutable std::vector<double> k3_;
  mutable std::vector<double> yt_;
};

} // namespace integrator

#endif // INTEGRATOR_RUNGE_KUTTA_3RD_HPP
