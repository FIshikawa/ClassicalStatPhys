#ifndef INTEGRATOR_POSITION_VELRET_HPP
#define INTEGRATOR_POSITION_VELRET_HPP

#include <string>
#include <vector>

namespace integrator {

class PositionVelret{
public:
  static std::string name() { return "Position Verlet method"; }
  PositionVelret(unsigned int dim) : dim_(dim), k1_(dim), k2_(dim), k3_(dim), yt_(dim) {}
  template<class F>
  void step(double t, double h, std::vector<double>& y, F const& f) const {
    const double h2 = h / 2;
    const int  dim2 = dim_ / 2;
    f(t, y , k1_);
    for (int i = 0; i < dim2; ++i) yt_[i] = y[i] + h2 * k1_[i];
    f(t + h, yt_, k2_);
    for (int i = dim2; i < dim_; ++i) y[i] = y[i] + h * k2_[i];
    f(t + h, y, k3_);
    for (int i = 0; i < dim2; ++i) y[i] = y[i] + h2 * k1_[i] + h2 * k3_[i];
  }
private:
  unsigned int dim_;
  mutable std::vector<double> k1_;
  mutable std::vector<double> k2_;
  mutable std::vector<double> k3_;
  mutable std::vector<double> yt_;
};

} // namespace integrator

#endif // INTEGRATOR_POSITION_VELRET_HPP
