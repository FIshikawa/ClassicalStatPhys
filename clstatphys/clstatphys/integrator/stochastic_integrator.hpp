#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include <string>
#include <vector>

namespace integrator {

class euler { //Euler Maruyama scheme
public:
  static std::string name() { return "Euler method"; }
  euler(double T, double gamma, unsigned int dim) : dim_(dim), k_(dim), T_(T), gamma_(gamma) {}
  template<class F>
  void step(double t, double h, std::vector<double>& y, F const& f, Rand & mt) const {
    boost::normal_distribution<> W(0.0,1.0);
    double D = std::sqrt(2 * gamma_ * T_ / dt);
    f(t, y, k_);
    for (int i = 0; i < dim_; ++i) y[i] = h * k_[i];
  }
private:
  unsigned int dim_;
  mutable std::vector<double> k_;
  double T_;
  double gamma_;
};//euler scheme end
  
class rk2 {//Heun scheme
public:
  static std::string name() { return "2nd-order Runge-Kutta method"; }
  rk2(unsigned int dim) : dim_(dim), k1_(dim), k2_(dim), yt_(dim) {}
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
};//rk2 end

} // namespace integrator

#endif // INTEGRATOR_HPP
