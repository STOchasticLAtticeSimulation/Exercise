#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// ------------ User may change --------------------------
const double sigma = 0.1; // coarse-graining param. 
const double kdx = 0.1; // ksigma Deltax / 2pi.
const int NL = 32;
const double Nf = 5.5;
const double dN = 0.01;
const int EXPSTEP = 10;
const std::string filename = "noisedata/largenoisetest.dat";
// -------------------------------------------------------

std::vector<double> dwlist(double N, double dN); 

// random distribution
std::random_device seed;
std::mt19937 engine(seed());
std::normal_distribution<> dist(0., 1.);

// useful macro
#define LOOP for(int i = 0; i < NL; i++) for(int j = 0; j < NL; j++) for(int k = 0; k < NL; k++)


int main()
{
  // ---------- start timer ----------
  struct timeval Nv;
  struct timezone Nz;
  double before, after;
  
  gettimeofday(&Nv, &Nz);
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  // --------------------------------------

#ifdef _OPENMP
  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;
#endif

  std::ofstream ofs(filename);
  if (ofs.fail()) {
    std::cout << "Data file coulnd't be opened." << std::endl;
    return -1;
  }
  
  std::vector<std::vector<double>> dwdata;
  double N = 0;
  int numsteps = 0;

  while (N<Nf+dN/2.) {
    if (numsteps % EXPSTEP == 0 && N < Nf+dN/2.) {
      std::cout << "\r" << std::setw(3) << N << "/" << Nf << std::flush;
    }
    dwdata.push_back(dwlist(N,dN));
    N += dN;
    numsteps++;
  }
  std::cout << std::endl;

  for (size_t i=0; i<dwdata[0].size(); i++) {
    for (size_t n=0; n<dwdata.size(); n++) {
      ofs << dwdata[n][i] << ' ';
    }
    ofs << std::endl;
  }
  

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  // -------------------------------------
}


std::vector<double> dwlist(double N,  double dN) {
  std::vector<double> dwlist(NL*NL*NL, 0);

  // make correlated noise
  double ksigma = 2. * M_PI * sigma * exp(N) / NL;
  double dtheta = 2. * M_PI * kdx / NL / ksigma;
  int divth = int(M_PI / dtheta);
  std::vector<double> divph(divth);
  std::vector<double> thetai(divth);
  std::vector<double> dphi(divth);
  std::vector<double> dOmegai(divth);

  // ------------ make noise ----------------
  std::vector<std::vector<double>> Omegalist;
  for (int n = 0; n < divth; n++) {
    thetai[n] = (n + 0.5) * dtheta;
    dphi[n] = 2. * M_PI * kdx / NL / ksigma / sin(thetai[n]);
    divph[n] = int(2. * M_PI / dphi[n]);
    for (int l = 0; l < divph[n]; l++){
      double dOmegai = sin(thetai[n]) * dtheta * dphi[n];
      double phii = l * dphi[n];
      Omegalist.push_back({thetai[n], dphi[n], phii, sqrt(dN) * dist(engine)});
    } 
  }
  // -----------------------------------------
  
  // -------------- add noise ----------------
#ifdef _OPENMP
#pragma omp parallel for
#endif
  LOOP{
    dwlist[i*NL*NL + j*NL + k] = 0; // initialize
    for (size_t n = 0; n < Omegalist.size(); n++){
      double thet = Omegalist[n][0];
      double dpht = Omegalist[n][1];
      double phit = Omegalist[n][2];
      double dwt = Omegalist[n][3];
      double dOmegai = sin(thet) * dtheta * dpht;
      double ksx = i * sin(thet)*cos(phit) + j * sin(thet)*sin(phit) + k * cos(thet);
      ksx *= ksigma;
      dwlist[i*NL*NL + j*NL + k] += 0.5 * sqrt(dOmegai / M_PI) * (cos(ksx) - sin(ksx)) * dwt;
    }
  }
  // -----------------------------------------
  
  return dwlist;
}
