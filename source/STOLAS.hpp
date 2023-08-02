#ifndef INCLUDED_STOLAS_
#define INCLUDED_STOLAS_

#include <sstream>
#include <fstream>
#include <iostream>
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
  const double Nprec = 0.01;
  // ----------------------------------------

  const std::string noisefilename = "noisetest.dat";

  std::string model;
  int NL;
  double Nf, dN, bias, Nbias, dNbias;
  std::ifstream noisefile;
  std::vector<std::vector<double>> noisedata;

public:
  STOLAS(){}
  STOLAS(std::string Model, double NF, std::string noisedir, int noisefileNo, double Bias, double NBias, double DNbias);

  double VV(double phi);
  double Vp(double phi);
  
  //void dNmap();

  double ep(double phi, double pi);
  double hubble(double phi, double pi);

  std::vector<double> dphidN(double N, std::vector<double> phi);
  std::vector<double> dphidNbias(double N, std::vector<double> phi, std::vector<double> pos);

  void RK4M(std::function<std::vector<double>(double, std::vector<double>)> dphidN, double dw, double &N, std::vector<double> &phi, double dN);
  void RK4Mbias(std::function<std::vector<double>(double, std::vector<double>, std::vector<double>)> dphidNbias, double dw, double &N, std::vector<double> &phi, double dN, std::vector<double> pos);
};

#endif
