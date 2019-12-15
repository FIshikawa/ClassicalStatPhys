#ifndef INTEGRATOR_VELOCITY_VELRET_PARALLEL_HPP
#define INTEGRATOR_VELOCITY_VELRET_PARALLEL_HPP

#include <string>
#include <vector>

namespace integrator {

class VelocityVelretParallel{
public:
  static std::string name() { return "Velocity Verlet method : parallelized"; }
  VelocityVelretParallel(unsigned int dim ) : dim_(dim), k1_(dim), k2_(dim) {}
  template<class F>
  void step(double t, double h, std::vector<double>& y, F const& f) const {
    const double h22 = h * h / 2;
    const double h2 = h / 2;
    const int  dim2 = dim_ / 2;
    const int dim3 = dim_ - dim2;
    f(t, y, k1_);
    #pragma omp parallel for 
    for (int i = 0; i < dim2; ++i){
      y[i] = y[i] +  h * k1_[i] + h22 * k1_[i + dim2];
    }
    f(t + h, y, k2_);
    #pragma omp parallel for 
    for (int i = 0; i < dim3; ++i){
      y[i+dim2] = y[i+dim2] +  h2 * (k2_[i+dim2]+k1_[i+dim2]);
    }
  }
private:
  unsigned int dim_;
  mutable std::vector<double> k1_;
  mutable std::vector<double> k2_;
};

} // namespace integrator

#endif // INTEGRATOR_VELOCITY_VELRET_PARALLEL_HPP
