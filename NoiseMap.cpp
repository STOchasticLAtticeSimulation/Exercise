// g++ -std=c++11 -O2 -o NoiseMap NoiseMap.cpp

#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sys/time.h>
#include <random>

using namespace std;

// Parameters
const int N = 17; // Number of lattice
const double L = (double) N;
double dt = 1.e-0;

// for random distribution
random_device seed;
mt19937 engine(seed());
normal_distribution<> dist(0., 1.);

// useful macro
#define LOOP for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) for(int k = 0; k < N; k++)


// Function prototype declaration
// vector<vector<vector<vector<vector<double>>>>> noise(double t, const vector<vector<vector<vector<vector<double>>>>> &x);


int main()
{
  // ---------------- start stop watch -----------------
  struct timeval tv;
  struct timezone tz;
  double before, after;

  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  // ---------------------------------------------------

  vector<vector<double>> xi{{0,0}}; // initial value
  vector<vector<vector<double>>> x1(N, xi); // dim 1
  vector<vector<vector<vector<double>>>> x2(N, x1); // dim 2
  vector<vector<vector<vector<vector<double>>>>> x(N, x2); // dim 3

  ofstream ofs_n("chaotic_noise.dat");

  vector<vector<vector<vector<vector<double>>>>> dxdw = x;

  // make correlation noise
  double peak = 4.; // peak point of power spectrum
  double ksigma = peak * 2. * M_PI / L; // * sigma * exp(t)
  double dtheta = M_PI / ksigma / L / 5.;
  int divth = int(M_PI / dtheta);
  vector<double> divph(divth);
  vector<double> thetai(divth);
  vector<double> dphi(divth);
  vector<double> dOmegai(divth);

  // inner product of coarse-grained vector and position vector
  double ksx = 0.;
  vector<vector<double>> Omegalist;
  for (int n = 0; n < divth; n++) {
    thetai[n] = (n + 0.5) * dtheta;
    dphi[n] = M_PI / ksigma / L / sin(thetai[n]) / 5.;
    divph[n] = int(2. * M_PI / dphi[n]);
    for (int l = 0; l < divph[n]; l++){
      double phii = l * dphi[n];
      Omegalist.push_back({thetai[n], dphi[n], phii, sqrt(dt) * dist(engine)});
    }
  }

  LOOP{
    dxdw[i][j][k][0][0] = 0.; // initialize
    for (int n = 0; n < Omegalist.size(); n++){
      double thet = Omegalist[n][0];
      double dpht = Omegalist[n][1];
      double phit = Omegalist[n][2];
      double dwt = Omegalist[n][3];
      double dOmegai = sin(thet) * dtheta * dpht;
      ksx = (i+0.5) * sin(thet)*cos(phit) + (j+0.5) * sin(thet)*sin(phit) + (k+0.5) * cos(thet);
      ksx *= ksigma;
      dxdw[i][j][k][0][0] += 0.5 * sqrt(dOmegai / M_PI) * (cos(ksx) - sin(ksx)) * dwt;
    }
  ofs_n << dxdw[i][j][k][0][0] << endl;
  }

  // ---------------- return elapsed time --------------
  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
  // ---------------------------------------------------
}


// vector<vector<vector<vector<vector<double>>>>> noise(double t, const vector<vector<vector<vector<vector<double>>>>> &x)
// {
  // vector<vector<vector<vector<vector<double>>>>> dxdw = x;

  // // make correlation noise
  // double ksigma = sigma * exp(t) * m;
  // double dtheta = M_PI / ksigma / L / 5.;
  // int divth = int(M_PI / dtheta);
  // vector<double> divph(divth);
  // vector<double> thetai(divth);
  // vector<double> dphi(divth);
  // vector<double> dOmegai(divth);
  // cout << divth << endl;

  // // inner product of coarse-grained vector and position vector
  // double ksx = 0.;
  // vector<vector<double>> Omegalist;
  // for (int n = 0; n < divth; n++) {
  //   thetai[n] = (n + 0.5) * dtheta;
  //   dphi[n] = M_PI / ksigma / L / sin(thetai[n]) / 5.;
  //   divph[n] = int(2. * M_PI / dphi[n]);
  //   for (int l = 0; l < divph[n]; l++){
  //     double phii = l * dphi[n];
  //     Omegalist.push_back({thetai[n], dphi[n], phii, sqrt(dt) * dist(engine)});
  //   }
  // }
  // cout << Omegalist.size() << endl;

  // LOOP{
  //   dxdw[i][j][k][0][0] = 0.; // initialize
  //   for (int n = 0; n < Omegalist.size(); n++){
  //     double thet = Omegalist[n][0];
  //     double dpht = Omegalist[n][1];
  //     double phit = Omegalist[n][2];
  //     double dwt = Omegalist[n][3];
  //     double dOmegai = sin(thet) * dtheta * dpht;
  //     ksx = (i+0.5) * sin(thet)*cos(phit) + (j+0.5) * sin(thet)*sin(phit) + (k+0.5) * cos(thet);
  //     ksx *= ksigma * dx;
  //     dxdw[i][j][k][0][0] += 0.5 * sqrt(dOmegai / M_PI) * (cos(ksx) - sin(ksx)) * dwt;
  //   }
  // cout << dxdw[i][j][k][0][0] << endl;
  // }

  // return dxdw;
// }
