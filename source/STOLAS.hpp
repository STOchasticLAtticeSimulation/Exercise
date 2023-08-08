#ifndef INCLUDED_STOLAS_
#define INCLUDED_STOLAS_

#define _USR_MATH_DEFINES
#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>

#ifdef _OPENMP
#include <omp.h>
#endif

class STOLAS
{
protected:
  // ------------ User may change -----------
  const double sigma = 0.1;
  const double kdx = 0.1;
  const double Nprec = 1e-4;
  // ----------------------------------------

  const std::string noisefilename = "noisetest.dat";
  const std::string Nfileprefix = "Nmap_Nf_";

  std::string model, Nfilename;
  int NL;
  double Nf, dN, bias, Nbias, dNbias;
  std::ifstream noisefile;
  std::vector<double> phii;
  std::vector<std::vector<double>> noisedata;

public:
  STOLAS(){}
  STOLAS(std::string Model, double NF, std::string noisedir, int noisefileNo, std::vector<double> Phii, double Bias, double NBias, double DNbias);

  double VV(double phi);
  double Vp(double phi);
  
  void dNmap();

  double ep(double phi, double pi);
  double hubble(double phi, double pi);

  std::vector<double> dphidN(double N, std::vector<double> phi);
  std::vector<double> dphidNbias(double N, std::vector<double> phi, int pos);

  void RK4(double &t, std::vector<double> &x, double dt);
  void RK4bias(double &t, std::vector<double> &x, double dt, int pos);
  void RK4M(double &N, std::vector<double> &phi, double dN, double dw);
  void RK4Mbias(double &N, std::vector<double> &phi, double dN, double dw, int pos);
};

#endif
