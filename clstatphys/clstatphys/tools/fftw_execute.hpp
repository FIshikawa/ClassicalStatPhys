#ifndef FFTW_EXECUTE_HPP
#define FFTW_EXECUTE_HPP

#include <iostream>
#include <complex>
#include <fftw3.h>

// todo : should modify structures 

namespace tools {


class FFTW { //todo : reveal its meaning : why did I implemet? 
public:
  SynmetrizeFFTW(int N): N_(N) {}
  void trans( std::vector<double>const& z, std::vector<double>& fz) const {
    fftw_complex *in, *out;
    fftw_plan p;
    in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*2*N_);
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*2*N_);
    p = fftw_plan_dft_1d(N_, in, out, FFTW_FORWARD,FFTW_ESTIMATE);
    
    for(int i = 0; i < N_; ++i){
      in[i][0] = z[i] ; 
      in[i+N_][0] = z[N_-i];
    }
    for(int i = 0;i < 2*N_; ++i){
      in[i][1] = 0.0;
    }
    fftw_execute(p);
    for(int i = 0; i < N_ ; ++i)fz[i] = out[i][0]/out[0][0] ;
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
  }
private:
  int N_;
}; //end ftrans

class PureFFTW{//non symetrize fourier trans
public:
  PureFFTW(int N): N_(N) {}
  void trans( std::vector<double>const& z, std::vector<double>& rfz, std::vector<double>& ifz ) const {
    fftw_complex *in, *out;
    fftw_plan p;
    in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*N_);
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*N_);
    p = fftw_plan_dft_1d(N_, in, out, FFTW_FORWARD,FFTW_ESTIMATE);
    
    for(int i = 0; i < N_; ++i){
      in[i][0] = z[i]; 
      in[i][1] = 0.0;
    }
    fftw_execute(p);
    for(int i = 0; i < N_; ++i){ 
    rfz[i] = out[i][0];
    ifz[i] = out[i][1];
    }
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
  }
private:
  int N_;
}; //end fftrans

class SynmetrizeFFTW { //synmetrize fourier trans. only real part
public:
  SynmetrizeFFTW(int N): N_(N) {}
  void trans( std::vector<double>const& z, std::vector<double>& fz) const {
    fftw_complex *in, *out;
    fftw_plan p;
    in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*N_);
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*N_);
    p = fftw_plan_dft_1d(N_, in, out, FFTW_FORWARD,FFTW_ESTIMATE);
    
    for(int i = 0; i < N_; ++i){
      in[i][0] = z[i]+z[N_+i] ; 
      in[i][1] = 0.0;
    }
    fftw_execute(p);
    for(int i = 0; i < N_ ; ++i)fz[i] = out[i][0]/out[0][0] ;
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
  }
private:
  int N_;
}; //end ftrans



}//end namespace

#endif //FFTW_EXECUTE_HPP
