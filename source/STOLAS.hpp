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
#include <complex>

#ifdef _OPENMP
#include <omp.h>
#endif

class STOLAS
{
protected:
  // ------------ User may change -----------
  const double sigma = 0.1;
  const double kdx = 0.1;
  const double Nprec = 1e-7;
  // ----------------------------------------

  const std::string noisefilename = "noisedata/noisemap_"; //"largenoisetest.dat";
  const std::string biasfilename = "biasdata/biasmap.dat";
  const std::string Nfileprefix = "data/Nmap_";
  const std::string Hfileprefix = "data/H_";
  const std::string pifileprefix = "data/pi_";
  const std::string wfileprefix = "data/logw_";
  const std::string powfileprefix = "data/power_";
  bool noisefilefail, biasfilefail;

  std::string model;
  int NL, noisefileNo;
  double dN, bias, Nbias, dNbias;
  std::ifstream noisefile, biasfile;
  std::ofstream Nfile, Hfile, pifile,wfile,powfile;
  std::vector<double> phii;
  std::vector<std::vector<double>> noisedata, biasdata, Hdata, pidata;
  std::vector<std::vector<std::vector<std::complex<double>>>> Nmap3D;

public:
  STOLAS(){}
  STOLAS(std::string Model, double DN, std::string sourcedir, int NoisefileNo, std::vector<double> Phii, double Bias, double NBias, double DNbias);

  bool checknoisefile();
  bool checkbiasfile();
  bool noisebiassize();
  bool Nfilefail();
  //bool Hfilefail();
  //bool pifilefail();
  bool wfilefail();
  bool powfilefail();
  
  double VV(double phi);
  double Vp(double phi);
  
  void dNmap();
  void animation();
  void powerspec();

  double ep(double phi, double pi);
  double hubble(double phi, double pi);

  std::vector<double> dphidN(double N, std::vector<double> phi);
  //std::vector<double> dphidNbias(double N, std::vector<double> phi, int pos);

  void RK4(double &t, std::vector<double> &x, double dt);
  //void RK4bias(double &t, std::vector<double> &x, double dt, int pos);
  //void RK4M(double &N, std::vector<double> &phi, double dN, double dw);
  void RK4Mbias(double &N, std::vector<double> &phi, double dN, double dw, double Bias //int pos
		);
};

#endif
