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
  // const std::string wfileprefix = "data/logw_";
  const std::string powfileprefix = "data/power_";
  const std::string cmpfileprefix = "data/compaction_";
  const std::string prbfileprefix = "data/probabilities";
  const std::string powsfileprefix = "data/powers";
  bool noisefilefail, biasfilefail;

  std::string model;
  int NL, noisefileNo;
  double dN, bias, Nbias, dNbias, phif;
  std::ifstream noisefile, biasfile;
  std::ofstream Nfile, Hfile, pifile, powfile, cmpfile, prbfile, powsfile; //, wfile
  std::vector<double> phii;
  std::vector<std::vector<double>> noisedata, biasdata, Hdata, pidata;
  std::vector<double> Ndata;
  std::vector<std::vector<std::vector<std::complex<double>>>> Nmap3D;

public:
  STOLAS(){}
  STOLAS(std::string Model, double DN, std::string sourcedir, int NoisefileNo, std::vector<double> Phii, double Bias, double NBias, double DNbias, double Phif);

  bool checknoisefile();
  bool checkbiasfile();
  bool noisebiassize();
  bool Nfilefail();
  //bool powfilefail();
  //bool cmpfilefail();
  
  double VV(double phi);
  double Vp(double phi);
  
  void dNmap();
  void animation();
  void powerspec();
  void compaction();

  double ep(double phi, double pi);
  double hubble(double phi, double pi);
  double calPphi(double &N, std::vector<double> &phi, double N0, bool broken);
  double calPpi(double &N, std::vector<double> &phi, double N0, bool broken);
  double RecalPphipi(double &N, std::vector<double> &phi, double N0, bool broken);
  bool EoI(std::vector<double> &phi);

  std::vector<double> dphidN(double N, std::vector<double> phi);

  void RK4(double &t, std::vector<double> &x, double dt);
  void RK4Mbias(double &N, std::vector<double> &phi, double dN, double dw, double Bias, double N0, bool broken);
};

#endif
