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
const std::string filename = "noisedata/noisemap_";
// -------------------------------------------------

const std::complex<double> II(0,1);

std::vector<double> dwlist(double N);
bool realpoint(int nx, int ny, int nz, int Num); // judge real point
bool complexpoint(int nx, int ny, int nz, int Num); // judge independent complex point
bool innsigma(int nx, int ny, int nz, int Num, double nsigma, double dn); // judge if point is in nsigma sphere shell

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
    std::cout << "The noise file couldn't be opened. 'mkdir noisedata'" << std::endl;
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
  std::vector<std::vector<double>> noisedata(totalstep, std::vector<double>(NL*NL*NL,0));
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<totalstep; i++) {
    noisedata[i] = dwlist(i*dN);

    /*#ifdef _OPENMP
#pragma omp critical
#endif
    {
      count++;
      std::cout << "\rNoiseGenerating : " << std::setw(3) << 100*count/totalstep << "%" << std::flush;
    }
    */
  }
  //  std::cout << std::endl;
  
  for (size_t i=0; i<noisedata[0].size(); i++) {
    for (size_t n=0; n<noisedata.size(); n++) {
      ofs << noisedata[n][i] << ' ';
    }
    ofs << std::endl;
    //std::cout << "\rExporting : " << std::setw(3) << 100*i/noisedata[0].size() << "%" << std::flush;
  }
  //  std::cout << "\rExporting : 100%" << std::endl;

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  // -------------------------------------
}


std::vector<double> dwlist(double N) {
  std::vector<std::vector<std::vector<std::complex<double>>>> dwk(NL, std::vector<std::vector<std::complex<double>>>(NL, std::vector<std::complex<double>>(NL, 0)));
  int count = 0;
  double nsigma = sigma*exp(N);
  std::vector<double> dwlist(NL*NL*NL,0);
  
  LOOP{
    if (innsigma(i,j,k,NL,nsigma,dn)) {
      if (realpoint(i,j,k,NL)) {
	dwk[i][j][k] = dist(engine);
	count++;
      } else if (complexpoint(i,j,k,NL)) {
	dwk[i][j][k] = (dist(engine) + II*dist(engine))/sqrt(2);
	count++;
      }
    }
  }

  // reflection
  int ip, jp, kp; // reflected index
  LOOP{
    if (innsigma(i,j,k,NL,nsigma,dn)) {
      if (!(realpoint(i,j,k,NL)||complexpoint(i,j,k,NL))) {
	if (i==0) {
	  ip = 0;
	} else {
	  ip = NL-i;
	}
	
	if (j==0) {
	  jp = 0;
	} else {
	  jp = NL-j;
	}
	
	if (k==0) {
	  kp = 0;
	} else {
	  kp = NL-k;
	}
	
	dwk[i][j][k] = conj(dwk[ip][jp][kp]);
	count++;
      }
    }
  }

  if (count==0) {
    return dwlist;
  }
  dwk /= sqrt(count);

  std::vector<std::vector<std::vector<std::complex<double>>>> dwlattice = fft(dwk);
  LOOP{
    dwlist[i*NL*NL + j*NL + k] = dwlattice[i][j][k].real();
  }

  return dwlist;
}

bool realpoint(int nx, int ny, int nz, int Num) {
  return (nx==0||nx==Num/2) && (ny==0||ny==Num/2) && (nz==0||nz==Num/2);
}

bool complexpoint(int nx, int ny, int nz, int Num) {
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

  return (1<=nxt && nxt!=Num/2 && nyt!=Num/2 && nzt!=Num/2) ||
    (nxt==Num/2 && nyt!=Num/2 && 1<=nzt && nzt!=Num/2) || (nxt!=Num/2 && 1<=nyt && nyt!=Num/2 && nzt==Num/2) || (1<=nxt && nxt!=Num/2 && nyt==Num/2 && nzt!=Num/2) ||
    (nxt==0 && nyt!=Num/2 && 1<=nzt && nzt!=Num/2) ||
    (nxt==Num/2 && nyt==Num/2 && 1<=nzt && nzt!=Num/2) || (nxt==Num/2 && 1<=nyt && nyt!=Num/2 && nzt==Num/2) || (1<=nxt && nxt!=Num/2 && nyt==Num/2 && nzt==Num/2) ||
    (nxt==0 && 1<=nyt && nyt!=Num/1 && nzt==0) ||
    (nxt==Num/2 && 1<=nyt && nyt!=Num/2 && nzt==0) || (1<=nxt && nxt!=Num/2 && nyt==0 && nzt==Num/2) || (nxt==0 && nyt==Num/2 && 1<=nzt && nzt!=Num/2);
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

  return std::abs(sqrt(nxt*nxt + nyt*nyt + nzt*nzt) - nsigma) <= dn/2.;
}

