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
