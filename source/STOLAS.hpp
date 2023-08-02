#ifndef INCLUDED_STOLAS_

#define INCLUDED_STOLAS_

#include <random>
#include "RK4.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

class STOLAS
{
protected:
  // ------------ User may change -----------
  const double sigma = 0.1;
  const double kdx = 0.1;
  // ----------------------------------------

  const std::string noisedataname = "noisetest.dat";

public:
  STOLAS(){}
  STOLAS(string model, double Nf, string noisedir, int noisedata, double bias, double Nbias, double DNbias);
  void dNmap();

  double VV(double phi);
  double Vp(double phi);
}
