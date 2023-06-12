// Inflection model　に importance sampling を実行．全ての omega に対して Normal Distributionを変更．
// g++ -std=c++11 -O2 -o inflection_IS_all inflection_IS_all.cpp
// 並列化するなら g++ -std=c++11 -O2 -Xpreprocessor -fopenmp -lomp -o inflection_IS_all inflection_IS_all.cpp
// g++ -std=c++11 -O2 -Xpreprocessor -fopenmp -L/opt/homebrew/opt/libomp/lib -lomp -o inflection_IS_all inflection_IS_all.cpp

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sys/time.h>
#include <functional>
#include <random>
#include <iomanip>
#include "MT.h"
#include "vec_op.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;


void RK4(function<vector<double>(double, vector<double>)> dphidN,double &N, vector<double> &phi, double dN); // e-folds N と (phi, pi) を渡すと EoM に従い dN だけ N, (phi, pi) を更新する
vector<double> dphidN(double N, vector<double> phi); // N と (phi, pi) の関数としての EoM
double ep(double phi, double pi);
double hubble(double phi, double pi);
double VV(double phi);
double Vp(double phi); // ポテンシャル VV の phi 微分
double Ncl(vector<double> phi,double N,double Nprec); // 初期条件 phi & N から end of inf までの e-folds Ncl を精度 Nprec で求める
 
// ------------ パラメータ ----------------- //
const string filename = "Ncl_inflection_biased.dat"; // 出力ファイル名(Ncl)
const string filename_c = "Lattice_inflection_biased.dat"; // 出力ファイル名(曲率ゆらぎ)
const string filename_f = "field_inflection_biased.dat"; // 出力ファイル名(phi, pi)
const double Nf = 5.0;  // lattice 終了時刻
const double dN = 0.01; // 時間刻み
const double AW = 0.02;
const double BW = 1;
const double CW = 0.04;
const double DW = 0;
const double GW = 3.076278e-2;
const double RW = 0.7071067;
const double CALV = 1000;
const double W0 = 12.35;//320182;//
const double CUP = 0.0382;
const double NPREC = 1e-7; // Ncl の精度
// ----------------------------------------- //

// 変数の初期値
// double N = 0.;// e-foldings
const double phi0 = 3.60547;
const double pi0 = -2.37409e-7;
const int NL = 17; // Number of lattice
const int N3 = NL * NL * NL; // for conveniensce
const double Ninv = 1. / NL; // for conveniensce
const double sigma = 1./10.; // coarse-grained scale parameter
const vector<double> xi{phi0, pi0}; // initial value
const double Nbias = 2.0; // Biased time

// output
//int numsteps = 0;

// random distribution
random_device seed;
mt19937 engine(seed());
normal_distribution<> dist(0., 1.);
normal_distribution<> dist1(5., 1.);

// useful macro
#define LOOP for(int i = 0; i < NL; i++) for(int j = 0; j < NL; j++) for(int k = 0; k < NL; k++)

// Euler-Maruyama method
template <class T>
void EulerM(function<T(double, const T&)> dphidNlist, function<T(double, const T&)> dwdNlist, double &N, T &x, double dN);

// for convenience
//vector<vector<double>> x1(NL, xi); // dim 1
//vector<vector<vector<double>>> x2(NL, x1); // dim 2
//vector<double> v1(NL,0);
//vector<vector<double>> v2(NL,v1);

vector<vector<vector<vector<double>>>> dwdNlist(double N, vector<vector<vector<vector<double>>>> xif);
vector<vector<vector<vector<double>>>> dphidNlist(double N, vector<vector<vector<vector<double>>>> xif);

double VV(double phi) {
  return W0*W0/CALV/CALV/CALV * ( CUP/pow(CALV,1./3) + AW/(exp(phi/sqrt(3))-BW) - CW/exp(phi/sqrt(3))
				  + exp(2*phi/sqrt(3))/CALV*(DW-GW/(RW*exp(sqrt(3)*phi)/CALV+1)) );
}

double Vp(double phi) {
  return exp(-phi/sqrt(3)) * ( CW - AW*exp(2*phi/sqrt(3))/(BW-exp(phi/sqrt(3)))/(BW-exp(phi/sqrt(3)))
			       + 3*exp(2*sqrt(3)*phi)*GW*RW/(CALV+exp(sqrt(3)*phi)*RW)/(CALV+exp(sqrt(3)*phi)*RW)
			       + 2*exp(sqrt(3)*phi)*(DW-CALV*GW/(CALV+exp(sqrt(3)*phi)*RW))/CALV )
    * W0*W0/sqrt(3)/CALV/CALV/CALV;
}


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
  cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << endl;
#endif

  // output
  int numsteps = 0;
  
  ofstream ofs(filename);
  ofstream ofs_c(filename_c);
  ofstream ofs_f(filename_f);
  ofs << setprecision(10);
  ofs_c << setprecision(10);
  ofs_f << setprecision(10);

  vector<vector<double>> x1(NL, xi); // dim 1
  vector<vector<vector<double>>> x2(NL, x1); // dim 2
  vector<double> v1(NL,0);
  vector<vector<double>> v2(NL,v1);
  vector<vector<vector<vector<double>>>> x(NL, x2); // dim 3

  // calculate average energy
  vector<vector<vector<double>>> de(NL, v2); // dim 3
  double average_e;

  double N = 0.;// e-foldings

  while (N < Nf) {
    // save the data
    if(numsteps % 10 == 0 && N < Nf){
      average_e = 0.; // initialize
      LOOP average_e += 3. * hubble(x[i][j][k][0], x[i][j][k][1]) * hubble(x[i][j][k][0], x[i][j][k][1]);
      average_e /= N3; // 平均エネルギー
      cout << N << ' ' << x[0][0][0][0] << ' ' << x[0][0][0][1] << ' ';
      vector<double> phi(2);
      LOOP{
        phi[0] = x[i][j][k][0];
        phi[1] = x[i][j][k][1];
        de[i][j][k] = 3. * hubble(phi[0], phi[1]) * hubble(phi[0], phi[1]) - average_e;
        ofs_c << de[i][j][k] / phi[1] / phi[1] / 3. << ' ';
      }
      cout << de[0][0][0]/phi[1]/phi[1]/3. << endl;
      ofs_c << endl;
    }
    numsteps++;

    EulerM<vector<vector<vector<vector<double>>>>>(dphidNlist, dwdNlist, N, x, dN); // Euler-Maruyama 1step <vector<vector<vector<vector<double>>>>>
  }

  LOOP ofs_f << x[i][j][k][0] << ' ' << x[i][j][k][1] << endl; // 最終的な場の値を出力
  

  //------------- Ncl -------------
  //*double N = Nf;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  LOOP{
    vector<double> phi{x[i][j][k][0], x[i][j][k][1]};
    
    double NCL = Ncl(phi,N,NPREC);
#ifdef _OPENMP
#pragma omp critical
#endif
    ofs // << setprecision(10)
      << i+NL*j+NL*NL*k << ' ' << NCL << endl;
  }

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
  // -------------------------------------
}

vector<vector<vector<vector<double>>>> dphidNlist(double N, vector<vector<vector<vector<double>>>> xif){
  vector<vector<vector<vector<double>>>> dphidNlist = xif; // Latticeと大きさを揃える
  LOOP{
    vector<double> phi(2);
    phi[0] = xif[i][j][k][0];
    phi[1] = xif[i][j][k][1];
    dphidNlist[i][j][k][0] = dphidN(N, phi)[0] * dN;
    dphidNlist[i][j][k][1] = dphidN(N, phi)[1] * dN;
  }

  return dphidNlist;
}

vector<vector<vector<vector<double>>>> dwdNlist(double N, vector<vector<vector<vector<double>>>> xif) {
  vector<vector<vector<vector<double>>>> dwdN = xif; // Latticeと大きさを揃える

  // make correlation noise
  double ksigma = 2. * M_PI * sigma * exp(N) * Ninv;
  double dtheta = 0.2 * M_PI * Ninv / ksigma;
  int divth = int(M_PI / dtheta);
  vector<double> divph(divth);
  vector<double> thetai(divth);
  vector<double> dphi(divth);
  vector<double> dOmegai(divth);

  // inner product of coarse-grained vector and position vector
  // double ksx = 0.;
  vector<vector<double>> Omegalist;
  if (Nbias <= N && N < Nbias+dN){
    cout << "biased!" << endl;
    for (int n = 0; n < divth; n++) {
      thetai[n] = (n + 0.5) * dtheta;
      dphi[n] = 0.2 * M_PI * Ninv / ksigma / sin(thetai[n]);
      divph[n] = int(2. * M_PI / dphi[n]);
      for (int l = 0; l < divph[n]; l++){
        double phii = l * dphi[n];
        Omegalist.push_back({thetai[n], dphi[n], phii, sqrt(dN) * dist1(engine)});
      }
    }
  } else {
    for (int n = 0; n < divth; n++) {
      thetai[n] = (n + 0.5) * dtheta;
      dphi[n] = 0.2 * M_PI * Ninv / ksigma / sin(thetai[n]);
      divph[n] = int(2. * M_PI / dphi[n]);
      for (int l = 0; l < divph[n]; l++){
        double phii = l * dphi[n];
        Omegalist.push_back({thetai[n], dphi[n], phii, sqrt(dN) * dist(engine)});
      }
    }
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  LOOP{
    dwdN[i][j][k][0] = 0.; // initialize
    dwdN[i][j][k][1] = 0.;
    for (int n = 0; n < Omegalist.size(); n++){
      double thet = Omegalist[n][0];
      double dpht = Omegalist[n][1];
      double phit = Omegalist[n][2];
      double dwt = Omegalist[n][3];
      double dOmegai = sin(thet) * dtheta * dpht;
      double ksx = (i+0.5) * sin(thet)*cos(phit) + (j+0.5) * sin(thet)*sin(phit) + (k+0.5) * cos(thet);
      ksx *= ksigma;
      dwdN[i][j][k][0] += 0.5 * sqrt(dOmegai / M_PI) * (cos(ksx) - sin(ksx)) * dwt;
    }
    dwdN[i][j][k][0] *= 0.5 * hubble(xif[i][j][k][0], xif[i][j][k][1]) / M_PI; // 係数を追加
    // if (Nbias <= N && N < Nbias+dN && i==0) cout << i << ' ' << j << ' ' << k << ' ' << hubble(xif[i][j][k][0], xif[i][j][k][1]) << ' ' << dwdN[i][j][k][0] << endl; // ノイズ確認用
  }

  return dwdN;
}

template <class T>
void EulerM(function<T(double, const T&)> dphidNlist, function<T(double, const T&)> dwdNlist, double &N, T &x, double dN)
{
  T xem = x;
  x += dphidNlist(N, xem);
  x += dwdNlist(N, xem);
  N += dN;
}

double Ncl(vector<double> phi,double N,double Nprec){
  // ofstream ofs(filename);
  double dN1 = dN;
  vector<double> prephi(2);

  while(dN1 >= Nprec) {
    while(ep(phi[0],phi[1])<=1.0){
      prephi[0]=phi[0];
      prephi[1]=phi[1];
      //ofs<< setprecision(10)<<N<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep(prephi[0],prephi[1])<<endl;
      RK4(dphidN, N, phi, dN1);
    // cout << setprecision(10) << ep(phi[0],phi[1]) << ' ' << phi[0] << ' ' << phi[1] << endl;
    }
    N -= dN1;
    
    phi[0]=prephi[0];
    phi[1]=prephi[1];
    dN1*=0.1;
  }
  return N; //Ncl;
}

double hubble(double phi, double pi) {
  return sqrt((pi*pi/2. + VV(phi))/3.);
}

double ep(double phi, double pi) {
  double HH = hubble(phi,pi);
  return pi*pi/2./HH/HH;
}

vector<double> dphidN(double N, vector<double> phi) {
  vector<double> dphidN(2);

  double xx = phi[0]; // phi
  double pp = phi[1]; // pi
  double HH = hubble(xx,pp);

  dphidN[0] = pp/HH;
  dphidN[1] = -3*pp - Vp(xx)/HH;

  return dphidN;
}

void RK4(function<vector<double>(double, vector<double>)> dphidN, double &N, vector<double> &phi, double dN) {
  vector<double> kx[4]; // 4段分の勾配ベクトル kx
  double a[4][4],b[4],c[4]; // ブッチャー係数

  // -------------- kx, a, b, c の初期化 --------------- //
  for(int i=0;i<=3;i++){
    kx[i].assign(phi.size(), 0.);
    
    for(int j=0;j<=3;j++){
      a[i][j]=0.;
    }
  }
  
  a[1][0]=1./2;  a[2][1]=1./2;  a[3][2]=1.;
  b[0]=1./6;     b[1]=1./3;     b[2]=1./3;    b[3]=1./6; 
  c[0]=0;        c[1]=1./2;     c[2]=1./2;    c[3]=1;
  // ----------------------------------------------------- //
  

  vector<double> X = phi; // i段目の位置用のベクトル
    
  for(int i=0;i<=3;i++){
    X = phi; // X の初期化
    
    for(int j=0;j<=3;j++){
      for (size_t xi = 0, size = phi.size(); xi < size; xi++) {
	X[xi] += dN * a[i][j] * kx[j][xi];
      }

      kx[i] = dphidN(N+c[i]*dN, X); // EoM
    }
  }

  for (size_t xi = 0, size = phi.size(); xi < size; xi++) {
    phi[xi] += dN*(b[0]*kx[0][xi] + b[1]*kx[1][xi] + b[2]*kx[2][xi] + b[3]*kx[3][xi]);
  }
  
  N+=dN;
}
