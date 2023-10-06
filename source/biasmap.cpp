#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <sys/time.h>

#include <chrono>
#include <thread>

#include "fft.hpp"
#include "vec_op.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


// -------------- User may change ------------------
const double sigma = 0.1; // ksigma = 2pi sigma exp(N) / L, nsigma = sigma exp(N)
const double dn = 1; // thickness of nsigma sphere shell
const int NL = pow(2,6); // box size L
const double dN = 0.01; // e-folds step
const std::string filename = "biasdata/biasmap.dat";
// -------------------------------------------------

const std::complex<double> II(0,1);

std::vector<double> biaslist(double N);
bool innsigma(int nx, int ny, int nz, int Num, double nsigma, double dn); // judge if point is in nsigma sphere shell

// random distribution
std::random_device seed;
std::mt19937 engine(seed());
std::normal_distribution<> dist(0., 1.);

// useful macro
#define LOOP for(int i = 0; i < NL; i++) for(int j = 0; j < NL; j++) for(int k = 0; k < NL; k++)


int main() 
{
  std::ofstream ofs(filename);
  if (ofs.fail()) {
    std::cout << "The bias file couldn't be opened." << std::endl;
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

  std::cout << "Box size : " << NL << std::endl;
  
  int totalstep = ceil(log((NL/2-1)/sigma)/dN), count = 0;
  std::vector<std::vector<double>> biasdata(totalstep, std::vector<double>(NL*NL*NL,0));
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<totalstep; i++) {
    biasdata[i] = biaslist(i*dN);
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      count++;
      std::cout << "\rBiasGenerating : " << std::setw(3) << 100*count/totalstep << "%" << std::flush;
    }
  }
  std::cout << std::endl;
  
  for (size_t i=0; i<biasdata[0].size(); i++) {
    for (size_t n=0; n<biasdata.size(); n++) {
      ofs << biasdata[n][i] << ' ';
    }
    ofs << std::endl;
    std::cout << "\rExporting : " << std::setw(3) << 100*i/biasdata[0].size() << "%" << std::flush;
  }
  std::cout << "\rExporting : 100%" << std::endl;

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  // -------------------------------------
}


std::vector<double> biaslist(double N) {
  std::vector<std::vector<std::vector<std::complex<double>>>> bk(NL, std::vector<std::vector<std::complex<double>>>(NL, std::vector<std::complex<double>>(NL, 0)));
  int count = 0;
  double nsigma = sigma*exp(N);
  std::vector<double> biaslist(NL*NL*NL,0);
  
  LOOP{
    if (innsigma(i,j,k,NL,nsigma,dn)) {
      bk[i][j][k] = 1;
      count++;
    }
  }

  if (count==0) {
    return biaslist;
  }
  bk /= count;

  std::vector<std::vector<std::vector<std::complex<double>>>> biaslattice = fft(bk);
  LOOP{
    biaslist[i*NL*NL + j*NL + k] = biaslattice[i][j][k].real();
  }

  return biaslist;
}

bool innsigma(int nx, int ny, int nz, int Num, double nsigma, double dn) {
  int nxt, nyt, nzt; // shifted index
  if (nx<=Num/2) {
    nxt = nx;
  } else {
    nxt = nx-Num;
  }

  if (ny<=Num/2) {
    nyt = ny;
  } else {
    nyt = ny-Num;
  }

  if (nz<=Num/2) {
    nzt = nz;
  } else {
    nzt = nz-Num;
  }

  return abs(sqrt(nxt*nxt + nyt*nyt + nzt*nzt) - nsigma) <= dn/2.;
}
