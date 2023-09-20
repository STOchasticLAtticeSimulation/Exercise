#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <sys/time.h>

#include <chrono>
#include <thread>

#include "fft.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


// -------------- User may change ------------------
const double sigma = 0.1; // ksigma = 2pi sigma exp(N) / L, nsigma = sigma exp(N)
const double dn = 1; // thickness of nsigma sphere shell
const int NL = pow(2,5); // box size L
const double dN = 0.01; // e-folds step
const std::string filename = "noisedata/noisemap_";
// -------------------------------------------------

std::vector<double> dwlist(double N);

// random distribution
std::random_device seed;
std::mt19937 engine(seed());
std::normal_distribution<> dist(0., 1.);

// useful macro
#define LOOP for(int i = 0; i < NL; i++) for(int j = 0; j < NL; j++) for(int k = 0; k < NL; k++)


int main(int argc, char* argv[]) 
{
  if (argc!=2) {
    std::cout << "Specify the noise file number correctly." << std::endl;
    return -1;
  }

  std::ofstream ofs(filename + std::string(argv[1]) + std::string(".dat"));
  if (ofs.fail()) {
    std::cout << "The noise file couldn't be opend." << std::endl;
    return -1;
  }

  
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

  int totalstep = ceil(log(NL/sigma)/dN), count = 0;
  std::vector<std::vector<double>> noisedata(totalstep, std::vector<double>(NL*NL*NL,0));
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<totalstep; i++) {
    noisedata[i] = dwlist(i*dN);
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      count++;
      std::cout << "\r" << std::setw(3) << 100*count/totalstep << "%" << std::flush;
    }
  }
  std::cout << std::endl;
}


std::vector<double> dwlist(double N) {
  std::this_thread::sleep_for(std::chrono::milliseconds(10));
  return std::vector<double>{1,2};
}

