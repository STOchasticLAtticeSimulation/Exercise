// g++ -std=c++11 -O2 -o NoiseMap NoiseMap.cpp

#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sys/time.h>
#include <random>
#include <functional>
#include "vec_op.hpp"

using namespace std;

// Parameters
#define m 1.e-5 // mass of phi
#define mm m*m // m^2
#define PHII 15.00 // initial value of phi
#define DPHII -0.1*m // initial value of phi

const int N = 17; // Number of lattice
const int N3 = N * N * N; // for conveniensce
const double Ninv = 1. / N; // for conveniensce
double t = 0.;
double tf = 5.5; // finish time
double dt = 1.e-2; // time interval
double sigma = 1./10.; // coarse-grained scale parameter
const string filename_c = "Lattice.dat"; // 出力ファイル名(曲率ゆらぎ)
const string filename_f = "field.dat"; // 出力ファイル名(phi, pi)
vector<double> xi{PHII, DPHII}; // initial value

// output
int numsteps = 0;

// random distribution
random_device seed;
mt19937 engine(seed());
normal_distribution<> dist(0., 1.);

// useful macro
#define LOOP for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) for(int k = 0; k < N; k++)

// Euler-Maruyama method
template <class T>
void EulerM(function<T(double, const T&)> dxdtlist, function<T(double, const T&)> dxdwlist, double &t, T &x, double dt);

// for convenience
vector<vector<double>> x1(N, xi); // dim 1
vector<vector<vector<double>>> x2(N, x1); // dim 2
vector<double> v1(N,0);
vector<vector<double>> v2(N,v1);

// Function prototype declaration
vector<double> dxdt(double t, const vector<double> phi);
double V(double phi);
double Vphi(double phi);
double H(double phi, double pi);
double epsilon(double phi, double pi);
vector<vector<vector<vector<double>>>> dxdwlist(double t, vector<vector<vector<vector<double>>>> xif);
vector<vector<vector<vector<double>>>> dxdtlist(double t, vector<vector<vector<vector<double>>>> xif);


int main()
{
  // ---------------- start stop watch -----------------
  struct timeval tv;
  struct timezone tz;
  double before, after;

  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  // ---------------------------------------------------

  ofstream ofs_c(filename_c);
  ofstream ofs_f(filename_f);
  ofs_c << setprecision(10);
  ofs_f << setprecision(10);
  vector<vector<vector<vector<double>>>> x(N, x2); // dim 3

  // calculate average energy
  vector<vector<vector<double>>> de(N, v2); // dim 3
  double average_e;

  while (t < tf) {
    average_e = 0.; // initialize
    LOOP average_e += 3. * H(x[i][j][k][0], x[i][j][k][1]) * H(x[i][j][k][0], x[i][j][k][1]);
    average_e /= N3; // average

    // save the data
    if(numsteps % 10 == 0 && t < tf){
      cout << t << ' ';
      vector<double> phi(2);
      LOOP{
        phi[0] = x[i][j][k][0];
        phi[1] = x[i][j][k][1];
        de[i][j][k] = 3. * H(phi[0], phi[1]) * H(phi[0], phi[1]) - average_e;
        ofs_c << de[i][j][k] / phi[1] / phi[1] / 3. << ' ';
      }
      cout << endl;
      ofs_c << endl;
    }
    numsteps++;

    EulerM<vector<vector<vector<vector<double>>>>>(dxdtlist, dxdwlist, t, x, dt); // Euler-Maruyama 1step
  }

  LOOP ofs_f << x[i][j][k][0] << ' ' << x[i][j][k][1] << endl; // 最終的な場の値を出力

  // ---------------- return elapsed time --------------
  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
  // ---------------------------------------------------
}

vector<double> dxdt(double t, vector<double> phi) {
  vector<double> dxdt(2); // initialize
  double xx = phi[0]; // phi
  double pp = phi[1]; // pi
  double HH = H(xx, pp);
  dxdt[0] = pp / HH; // definition of phi0dot
  dxdt[1] = - 3. * pp - Vphi(xx) / HH; // EoM
  dxdt *= dt;

  return dxdt;
}

double V(double phi) {
  return 0.5 * mm * phi * phi;
}

double Vphi(double phi) {
  return mm * phi;
}

double H(double phi, double pi) {
  return sqrt((0.5 * pi * pi + V(phi)) / 3.);
}

double epsilon(double phi, double pi) {
  return 0.5 * (pi * pi) / H(phi, pi) / H(phi, pi);
}

vector<vector<vector<vector<double>>>> dxdtlist(double t, vector<vector<vector<vector<double>>>> xif){
  vector<vector<vector<vector<double>>>> dxdtlist = xif; // Latticeと大きさを揃える
  LOOP{
    vector<double> phi(2);
    phi[0] = xif[i][j][k][0];
    phi[1] = xif[i][j][k][1];
    dxdtlist[i][j][k][0] = dxdt(t, phi)[0];
    dxdtlist[i][j][k][1] = dxdt(t, phi)[1];
  }

  return dxdtlist;
}

vector<vector<vector<vector<double>>>> dxdwlist(double t, vector<vector<vector<vector<double>>>> xif) {
  vector<vector<vector<vector<double>>>> dxdw = xif; // Latticeと大きさを揃える

  // make correlation noise
  double ksigma = 2. * M_PI * sigma * exp(t) * Ninv;
  double dtheta = 0.2 * M_PI * Ninv / ksigma;
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
    dphi[n] = 0.2 * M_PI * Ninv / ksigma / sin(thetai[n]);
    divph[n] = int(2. * M_PI / dphi[n]);
    for (int l = 0; l < divph[n]; l++){
      double phii = l * dphi[n];
      Omegalist.push_back({thetai[n], dphi[n], phii, sqrt(dt) * dist(engine)});
    }
  }

  LOOP{
    dxdw[i][j][k][0] = 0.; // initialize
    dxdw[i][j][k][1] = 0.;
    for (int n = 0; n < Omegalist.size(); n++){
      double thet = Omegalist[n][0];
      double dpht = Omegalist[n][1];
      double phit = Omegalist[n][2];
      double dwt = Omegalist[n][3];
      double dOmegai = sin(thet) * dtheta * dpht;
      ksx = (i+0.5) * sin(thet)*cos(phit) + (j+0.5) * sin(thet)*sin(phit) + (k+0.5) * cos(thet);
      ksx *= ksigma;
      dxdw[i][j][k][0] += 0.5 * sqrt(dOmegai / M_PI) * (cos(ksx) - sin(ksx)) * dwt;
    }
    dxdw[i][j][k][0] *= 0.5 * H(xif[i][j][k][0], xif[i][j][k][1]) / M_PI; // 係数を追加
  }

  return dxdw;
}


template <class T>
void EulerM(function<T(double, const T&)> dxdtlist, function<T(double, const T&)> dxdwlist, double &t, T &x, double dt)
{
  x += dxdtlist(t, x);
  x += dxdwlist(t, x);// 0.5 * H(x[0], x[1]) / M_PI * dxdw(t);
  t += dt;
}
