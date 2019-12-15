#ifndef INTEGRATOR_RUNGE_KUTTA_4TH_HPP
#define INTEGRATOR_RUNGE_KUTTA_4TH_HPP

#include <string>
#include <vector>

namespace integrator {

class RungeKutta4th{
public:
  static std::string name() { return "4th-order Runge-Kutta method"; }
  RungeKutta4th(int dim) : dim_(dim), k1_(dim), k2_(dim), k3_(dim), k4_(dim), yt_(dim) {}
  template<class F>
  void step(double t, double h, std::vector<double>& y, F const& f) const {
    double h2 = h / 2;
    f(t, y, k1_);
    for (int i = 0; i < dim_; ++i) yt_[i] = y[i] + h2 * k1_[i];
    f(t + h2, yt_, k2_);
    for (int i = 0; i < dim_; ++i) yt_[i] = y[i] + h2 * k2_[i];
    f(t + h2, yt_, k3_);
    for (int i = 0; i < dim_; ++i) yt_[i] = y[i] + h * k3_[i];
    f(t + h, yt_, k4_);
    for (int i = 0; i < dim_; ++i)
      y[i] = y[i] + (h / 6) * (k1_[i] + 2 * k2_[i] + 2 * k3_[i] + k4_[i]);
  }
private:
  int dim_;
  mutable std::vector<double> k1_;
  mutable std::vector<double> k2_;
  mutable std::vector<double> k3_;
  mutable std::vector<double> k4_;
  mutable std::vector<double> yt_;
};

} // namespace integrator

#endif // INTEGRATOR_RUNGE_KUTTA_4TH_HPP
