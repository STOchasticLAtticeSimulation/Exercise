// NoiseMap.cpp と Ncl.cpp を組み合わせ、end of inflation での N のマップ (or δN のマップ) を出力する lattice simulation code を作成。

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sys/time.h>
#include <functional>
#include <random>
#include <iomanip>
#include <string>
#include <functional>
#include "vec_op.hpp"
#include "MT.h"

using namespace std;

void RK4(function<vector<double>(double, vector<double>)> dphidN,double &N, vector<double> &phi, double dN); // e-folds N と (phi, pi) を渡すと EoM に従い dN だけ N, (phi, pi) を更新する
vector<double> dphidN(double N, vector<double> phi); // N と (phi, pi) の関数としての EoM
double ep(double phi, double pi);
double hubble(double phi, double pi);
double VV(double phi);
double Vp(double phi); // ポテンシャル VV の phi 微分
double Ncl(vector<double> phi,double N,double Nprec); // 初期条件 phi & N から end of inf までの e-folds Ncl を精度 Nprec で求める
vector<vector<vector<vector<vector<double>>>>> dxdw(double t, const vector<vector<vector<vector<vector<double>>>>> &x);//相関を持ったガウシアンノイズ

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
void EulerM(function<T(double, const T&)> dphidN, function<T(double, const T&)> dxdw, double &t, T &x, double dt);



// ------------ パラメータ ----------------- //
#define m 1.e-5 // mass of phi
#define PHII 15.00 // initial value of phi
const string filename = "lattice_simulation.dat"; // 出力ファイル名
const int NL = 17; // 一辺の格子数
const int N3 = NL * NL * NL; // 全格子数
const double Ninv = 1. / NL; // for conveniensce
const double N = 0.;// lattice 開始時刻
const double Nf = 5.5; // lattice 終了時刻
const double dN = 0.01; // 時間刻み
const double mm=1.0e-5; // 質量
const double NPREC = 1e-7; // Ncl の精度
// ----------------------------------------- //

// 変数の初期値を代入
double sigma = 1./10.; // coarse-grained scale parameter
vector<vector<double>> xi{{15.0,0.0}}; // 各点の phi と pi の初期値
vector<vector<vector<double>>> x1(NL, xi); // dim 1
vector<vector<vector<vector<double>>>> x2(NL, x1); // dim 2
vector<vector<vector<vector<vector<double>>>>> x(NL, x2); // dim 3

const double phi0;//　lattice 計算後のphiの値
const double pi0;// lattice 計算後のpiの値

double VV(double phi) {
  return mm*mm*phi*phi/2.;
}

double Vp(double phi) {
  return mm*mm*phi;
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

// lattice 計算始め
  // for calculate average energy
  // vector<vector<vector<double>>> de = hubble(&x);
  // double average_e;

  while (N < Nf) {
    // average_e = 0.; // initialize
    // LOOP average_e += 3. * hubble(x)[i][j][k] * hubble(x)[i][j][k];
    // average_e /= N3; // average

    // // save the data
    // if(numsteps % 10 == 0 && N < Nf){
    //   cout << N << ' ';
    //   LOOP{
    //     de[i][j][k] = 3. * hubble(x)[i][j][k] * hubble(x)[i][j][k] - average_e;
    //     //ofs_c << de[i][j][k] / x[i][j][k][0][1] / x[i][j][k][0][1] / 3. << ' ';
    //   }
    //   cout << endl;
    //   //ofs_c << endl;
    // }
    // numsteps++;

    EulerM<vector<vector<vector<vector<vector<double>>>>>>(dphidN, dwdN, N, x, dN); // Euler-Maruyama 1step
  }

  LOOP cout << x[i][j][k][0][0] << ' ' << x[i][j][k][0][1] << endl; // 最終的な場の値を出力



// Ncl 計算始め
  double N = Nf;
  vector<double> phi{phi0, pi0};

  double NCL=Ncl(phi,N,NPREC);
  cout<< setprecision(10)<<NCL<<endl;

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
  // -------------------------------------
}

double hubble(double phi, double pi) {
  return sqrt((pi*pi/2. + VV(phi))/3.);
}

double ep(double phi, double pi) {
  double HH = hubble(phi,pi);
  return pi*pi/2./HH/HH;
}

vector<vector<vector<vector<vector<double>>>>> dwdN(double t, const vector<vector<vector<vector<vector<double>>>>> &x)//相関を持ったガウシアンノイズ
{
  vector<vector<vector<vector<vector<double>>>>> dwdN = x;
  //vector<vector<vector<double>>> Hc = hubble(x[i][j][k][0][0],x[i][j][k][0][1]);

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
      Omegalist.push_back({thetai[n], dphi[n], phii, sqrt(dN) * dist(engine)});
    }
  }

  LOOP{
    dwdN[i][j][k][0][0] = 0.; // initialize
    dwdN[i][j][k][0][1] = 0.;
    for (int n = 0; n < Omegalist.size(); n++){
      double thet = Omegalist[n][0];
      double dpht = Omegalist[n][1];
      double phit = Omegalist[n][2];
      double dwt = Omegalist[n][3];
      double dOmegai = sin(thet) * dtheta * dpht;
      ksx = (i+0.5) * sin(thet)*cos(phit) + (j+0.5) * sin(thet)*sin(phit) + (k+0.5) * cos(thet);
      ksx *= ksigma;
      dwdN[i][j][k][0][0] += 0.5 * sqrt(dOmegai / M_PI) * (cos(ksx) - sin(ksx)) * dwt;
    }
    dwdN[i][j][k][0][0] *= 0.5 * hubble(x[i][j][k][0][0],x[i][j][k][0][1]) / M_PI;//ノイズ項
  }

  return dwdN;
}


void EulerM(function<vector<double>(double, vector<double>)> dphidN, function<vector<vector<vector<vector<vector<double>>>>>(double, const vector<vector<vector<vector<vector<double>>>>>)> dwdN, double &N,vector<vector<vector<vector<vector<double>>>>> &x, double dN)
{
  LOOP {x[i][j][k][0][0] += dphidN[0];
        x[i][j][k][0][1] += dphidN[1];
  }
  x += dwdN(N, x);
  N += dN;
}


vector<double> dphidN(double N, vector<double> phi) {
  vector<double> dphidN(2);

  double xx = phi[0]; // phi
  double pp = phi[1]; // pi
  double HH = hubble(xx,pp);

  //dphidN[0] = (phi[1]/(sqrt((1.0/6.0)*((pow(phi[1],2.0))+pow((mm*phi[0]),2.0)))));//dphi/dN=(pi/H)　ノイズなし
  //dphidN[1] = -3.0*phi[1]-((pow(mm,2.0))*phi[0])/(sqrt((1.0/6.0)*((pow(phi[1],2.0))+pow((mm*phi[0]),2.0))));//dpi/dN=-3pi-(v'/H)

  dphidN[0] = pp/HH;
  dphidN[1] = -3*pp - Vp(xx)/HH;

  return dphidN;
}

double Ncl(vector<double> phi,double N,double Nprec){
  ofstream ofs(filename);
  double dN1 = dN;
  vector<double> prephi(2);

  while(dN1 >= Nprec) {
    while(ep(phi[0],phi[1])<=1.0){
      prephi[0]=phi[0];
      prephi[1]=phi[1];
      ofs<< setprecision(10)<<N<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep(prephi[0],prephi[1])<<endl;
      RK4(dphidN, N, phi, dN1);
    }
    N -= dN1;
    
    phi[0]=prephi[0];
    phi[1]=prephi[1];
    dN1*=0.1;
  }
  return N; //Ncl;
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
