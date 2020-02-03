#ifndef STAT_ACCUMULATOR_HPP
#define STAT_ACCUMULATOR_HPP

#include <cmath>
#include <random>

namespace tools {

class Accumulator {
public:
  Accumulator() { reset(); }
  void reset() {
    count_ = 0;
    sum1_ = sum2_ = sum3_ = sum4_ = 0;
  }
  double operator<<(double v) {
    ++count_;
    sum1_ += v;
    sum2_ += v * v;
    sum3_ += v * v * v;
    sum4_ += v * v * v * v;
    return v;
  }

  Accumulator& operator+(const Accumulator & accumulator){
    this->sum1_ += accumulator.sum1();
    this->sum2_ += accumulator.sum2();
    this->sum3_ += accumulator.sum3();
    this->sum4_ += accumulator.sum4();
    this->count_ += accumulator.count();
    return *this;
  }

  void add(double sum1, double sum2, double sum3, double sum4, long count){
    sum1_ += sum1;
    sum2_ += sum2;
    sum3_ += sum3;
    sum4_ += sum4;
    count_ += count;
  }

  long count() const { return count_; }

  double sum1() const { return sum1_;}
  double sum2() const { return sum2_;}
  double sum3() const { return sum3_;}
  double sum4() const { return sum4_;}
  double moment1() const { return count() ? (sum1_ / count()) : 0; }
  double moment2() const { return count() ? (sum2_ / count()) : 0; }
  double moment3() const { return count() ? (sum3_ / count()) : 0; }
  double moment4() const { return count() ? (sum4_ / count()) : 0; }

  double central_moment1() const { return 0; }
  double central_moment2() const {
    double n = count();
    return (n > 1) ? (n/(n-1)) * (moment2() - moment1() * moment1()): 0;
  }
  double central_moment3() const {
    double n = count();
    double central_moment3_temp = moment3() - 3 * moment2() * moment1() + 2 * moment1() * moment1() * moment1();
    return (n > 2) ? (n*n/((n-1)*(n-2))) * central_moment3_temp : 0;
  }
  double central_moment4() const {
    double n = count();
    double central_moment4_temp = moment4() - 4 * moment3() * moment1() 
      + 6 * moment2() * moment1() * moment1() - 3 * moment1() * moment1() * moment1() * moment1();
    double square_variance_temp = (moment2() * moment2() - 2 * moment1() * moment1() * moment2() + moment1() * moment1() * moment1() * moment1());
    return (n > 3) ?
      n/((n-1)*(n-2)*(n-3)) *
      ((n*n-2*n+3)*central_moment4_temp - (6*n-9)*square_variance_temp) : 0;
  }
  double cumulant1() const { return moment1(); }
  double cumulant2() const { return central_moment2(); }
  double cumulant3() const { return central_moment3(); }
  double cumulant4() const {
    double n = count();
    double central_moment4_temp = moment4() - 4 * moment3() * moment1() 
      + 6 * moment2() * moment1() * moment1() - 3 * moment1() * moment1() * moment1() * moment1();
    double square_variance_temp = (moment2() * moment2() - 2 * moment1() * moment1() * moment2() + moment1() * moment1() * moment1() * moment1());
    return (n > 3) ?
      n*n/((n-1)*(n-2)*(n-3)) *
      ((n+1)*central_moment4_temp - 3*(n-1)*square_variance_temp) : 0;
  }
  double mean() const { return cumulant1(); }
  double average() const { return mean(); }
  double variance() const { return cumulant2(); }
  double square_variance() const {
    double n = count();
    double central_moment4_temp = moment4() - 4 * moment3() * moment1() 
      + 6 * moment2() * moment1() * moment1() - 3 * moment1() * moment1() * moment1() * moment1();
    double square_variance_temp = (moment2() * moment2() - 2 * moment1() * moment1() * moment2() + moment1() * moment1() * moment1() * moment1());
    return (n > 3) ?
      n/((n-1)*(n-2)*(n-3)) *
      ((n*n-3*n+3)*square_variance_temp-(n-1)*central_moment4_temp) : 0;
  }
  double variance_error() const {return (count_ > 3) ? std::sqrt(central_moment4()/count() - (count() - 3)/(count()*(count()-1)) * square_variance()) : 0; }
  double error() const { return (count_ > 1) ? std::sqrt(variance() / count()) : 0; }
  double standard_deviation() const { return std::sqrt(variance()); }
  double skewness() const {
    return (standard_deviation() > 0) ? central_moment3() / (standard_deviation() * standard_deviation() * standard_deviation()) : 0;
  }
  double kurtosis() const {
    return (standard_deviation() > 0) ? central_moment4() / (standard_deviation() * standard_deviation() * standard_deviation() * standard_deviation()) : 0;
  }
  double kurtosis_excess() const { return kurtosis() - 3; }

private:
  long count_;
  double sum1_, sum2_, sum3_, sum4_;
};

} // end namespace stat

#endif // STAT_ACCUMULATOR_HPP
