// g++ -std=c++11 -O2 -o NoiseMap NoiseMap.cpp

#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sys/time.h>
#include <random>
#include <functional>
#include <iomanip>
#include "vec_op.hpp"

using namespace std;

// Parameters
#define m 1.e-5 // mass of phi
#define mm m*m // m^2
#define PHII 15.00 // initial value of phi

const int N = 17; // Number of lattice
const int N3 = N * N * N; // for conveniensce
const double Ninv = 1. / N; // for conveniensce
double t = 0.;
double tf = 5.5; // finish time
double dt = 1.e-2; // time interval
double sigma = 1./10.; // coarse-grained scale parameter
const string filename_c = "Lattice.dat"; // 出力ファイル名(曲率ゆらぎ)
const string filename_f = "field.dat"; // 出力ファイル名(phi, pi)

// for output
int numsteps = 0;

// for random distribution
random_device seed;
mt19937 engine(seed());
normal_distribution<> dist(0., 1.);

// useful macro
#define LOOP for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) for(int k = 0; k < N; k++)

// Euler-Maruyama method
template <class T>
void EulerM(function<T(double, const T&)> dxdt, function<T(double, const T&)> dxdw, double &t, T &x, double dt);

// for convenience
vector<double> v1(N,0);
vector<vector<double>> v2(N,v1);

// Function prototype declaration
vector<vector<vector<vector<vector<double>>>>> dxdt(double t, const vector<vector<vector<vector<vector<double>>>>> &x);
vector<vector<vector<double>>> V(const vector<vector<vector<vector<vector<double>>>>> &x);
vector<vector<vector<double>>> Vphi(const vector<vector<vector<vector<vector<double>>>>> &x);
vector<vector<vector<double>>> H(const vector<vector<vector<vector<vector<double>>>>> &x);
vector<vector<vector<double>>> epsilon(const vector<vector<vector<vector<vector<double>>>>> &x);
vector<vector<vector<vector<vector<double>>>>> dxdw(double t, const vector<vector<vector<vector<vector<double>>>>> &x);


int main()
{
  // ---------------- start stop watch -----------------
  struct timeval tv;
  struct timezone tz;
  double before, after;

  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  // ---------------------------------------------------

  vector<vector<double>> xi{{15.0,-0.1*m}}; // initial value
  vector<vector<vector<double>>> x1(N, xi); // dim 1
  vector<vector<vector<vector<double>>>> x2(N, x1); // dim 2
  vector<vector<vector<vector<vector<double>>>>> x(N, x2); // dim 3
  ofstream ofs_c(filename_c);
  ofstream ofs_f(filename_f);
  ofs_c << setprecision(10);
  ofs_f << setprecision(10);

  // for calculate average energy
  vector<vector<vector<double>>> de = H(x);
  double average_e;

  while (t < tf) {
    average_e = 0.; // initialize
    LOOP average_e += 3. * H(x)[i][j][k] * H(x)[i][j][k];
    average_e /= N3; // average

    // save the data
    if(numsteps % 10 == 0 && t < tf){
      cout << t << ' ';
      LOOP{
        de[i][j][k] = 3. * H(x)[i][j][k] * H(x)[i][j][k] - average_e;
        ofs_c << de[i][j][k] / x[i][j][k][0][1] / x[i][j][k][0][1] / 3. << ' ';
      }
      cout << endl;
      ofs_c << endl;
    }
    numsteps++;

    EulerM<vector<vector<vector<vector<vector<double>>>>>>(dxdt, dxdw, t, x, dt); // Euler-Maruyama 1step
  }

  LOOP ofs_f << x[i][j][k][0][0] << ' ' << x[i][j][k][0][1] << endl; // 最終的な場の値を出力

  // ---------------- return elapsed time --------------
  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
  // ---------------------------------------------------
}

vector<vector<vector<vector<vector<double>>>>> dxdt(double t, const vector<vector<vector<vector<vector<double>>>>> &x)
{
  vector<vector<vector<vector<vector<double>>>>> dxdt = x; // initialize
  vector<vector<vector<double>>> Hc = H(x);
  vector<vector<vector<double>>> Vp = Vphi(x);
  
  LOOP{
    double phi0dot = x[i][j][k][0][1];
    dxdt[i][j][k][0][0] = phi0dot / Hc[i][j][k]; // definition of phi0dot
    dxdt[i][j][k][0][1] = - 3. * phi0dot - Vp[i][j][k] / Hc[i][j][k]; // EoM
  }
  dxdt *= dt;

  return dxdt;
}

vector<vector<vector<double>>> V(const vector<vector<vector<vector<vector<double>>>>> &x)
{
  vector<vector<vector<double>>>  v(N,v2);
  LOOP{
    double phi0 = x[i][j][k][0][0];
    v[i][j][k] = 0.5 * mm * phi0 * phi0;
  }

  return v;
}

vector<vector<vector<double>>> Vphi(const vector<vector<vector<vector<vector<double>>>>> &x)
{
  vector<vector<vector<double>>>  v(N,v2);
  LOOP{
    double phi0 = x[i][j][k][0][0];
    v[i][j][k] = mm * phi0;
  }

  return v;
}

vector<vector<vector<double>>> H(const vector<vector<vector<vector<vector<double>>>>> &x)
{
  vector<vector<vector<double>>>  v(N,v2);
  vector<vector<vector<double>>> Vc = V(x);
  LOOP{
    double phi0dot = x[i][j][k][0][1];
    v[i][j][k] = sqrt((phi0dot * phi0dot / 2. + Vc[i][j][k]) / 3. );
  }

  return v;
}

vector<vector<vector<double>>> epsilon(const vector<vector<vector<vector<vector<double>>>>> &x)
{
  vector<vector<vector<double>>>  v(N,v2);
  vector<vector<vector<double>>> Hc = H(x);
  LOOP{
    double phi0dot = x[i][j][k][0][1];
    v[i][j][k] = (phi0dot * phi0dot) / 2. / Hc[i][j][k] / Hc[i][j][k];
  }

  return v;
}

vector<vector<vector<vector<vector<double>>>>> dxdw(double t, const vector<vector<vector<vector<vector<double>>>>> &x)
{
  vector<vector<vector<vector<vector<double>>>>> dxdw = x;
  vector<vector<vector<double>>> Hc = H(x);

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
    dxdw[i][j][k][0][0] = 0.; // initialize
    dxdw[i][j][k][0][1] = 0.;
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
    dxdw[i][j][k][0][0] *= 0.5 * Hc[i][j][k] / M_PI;
  }

  return dxdw;
}


template <class T>
void EulerM(function<T(double, const T&)> dxdt, function<T(double, const T&)> dxdw, double &t, T &x, double dt)
{
  x += dxdt(t, x);
  x += dxdw(t, x); 
  t += dt;

  /* ------------------
    1行目で更新した x を2行目の dxdw に使ってしまっているのはまずい...?
    オイラー丸山なら OK?
  ------------------- */
}
