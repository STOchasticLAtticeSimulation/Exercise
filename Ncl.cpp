// oscillation.cpp を参考に chaotic potential に対して指定された初期条件 (phi, pi) から end of inflation までの e-folds を小数点以下 7 桁まで求めるコード Ncl.cpp を作成
//Nclを"Ncl.dat"に出力する

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sys/time.h>
#include <functional>
#include <random>
#include <iomanip>
#include "MT.h"

using namespace std;

void RK4(function<vector<double>(double, vector<double>)> pi,double &N, vector<double> &phi, double dN); // e-folds N と (phi, pi) を渡すと EoM に従い dN だけ N, (phi, pi) を更新する
vector<double> dphidN(double N, vector<double> phi); // N と (phi, pi) の関数としての EoM
vector<double> ep(vector<double> phi,vector<double> pi);
vector<double> hubble(vector<double> phi,vector<double> pi);
vector<double> VV(vector<double> phi);
vector<double> Vp(vector<double> phi); // ポテンシャル VV の phi 微分
 
// ------------ パラメータ ----------------- //
const string filename = "Ncl.dat"; // 出力ファイル名
const double Nf = 7.0; // lattice 終了時刻
const double dN = 0.01; // 時間刻み
const double mm=1.0e-5; // 質量
// ----------------------------------------- //

// 変数の初期値
const double phi0 = 14;
const double pi0 = -8.5e-6;

// 乱数
// double Uniform( void ){
// 	return genrand_real3();
// }

// double rand_normal(double mu,double sigma){
//   double z=sqrt(-2.0*log(Uniform()))*sin(2.0*M_PI*Uniform());
//   return mu+sigma*z;//sigmaは標準偏差。分散0.01にしたいときはsigma=0.1を入れる
// }

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
  
  
  ofstream ofs(filename);
  double N = Nf;
  vector<double> phi{phi0, pi0};
  vector<double> prephi(2);
  double ep=pow(phi[1],2.0)/2.0/(pow((sqrt((1.0/6.0)*((pow(phi[1],2.0))+pow((mm*phi[0]),2.0)))),2.0));
  
  while(ep<=1.0){
    prephi[0]=phi[0];
    prephi[1]=phi[1];
    ofs<< setprecision(10)<<N<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep<<endl;
    RK4(dphidN, N, phi, dN);
    ep=pow(phi[1],2.0)/2.0/(pow((sqrt((1.0/6.0)*((pow(phi[1],2.0))+pow((mm*phi[0]),2.0)))),2.0));
  }
  double Ncl=N-dN;
  ep=pow(prephi[1],2.0)/2.0/(pow((sqrt((1.0/6.0)*((pow(prephi[1],2.0))+pow((mm*prephi[0]),2.0)))),2.0));

  cout<< setprecision(10)<<Ncl<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep<<endl;
  
  phi[0]=prephi[0];
  phi[1]=prephi[1];
  double dN1=dN*0.1;
  N=Ncl;
  while(ep<=1.0){
    prephi[0]=phi[0];
    prephi[1]=phi[1];
    ofs<< setprecision(10)<<N<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep<<endl;
    RK4(dphidN, N, phi, dN1);
    ep=pow(phi[1],2.0)/2.0/(pow((sqrt((1.0/6.0)*((pow(phi[1],2.0))+pow((mm*phi[0]),2.0)))),2.0));
  }
  Ncl=N-dN1;
  ep=pow(prephi[1],2.0)/2.0/(pow((sqrt((1.0/6.0)*((pow(prephi[1],2.0))+pow((mm*prephi[0]),2.0)))),2.0));

  cout<< setprecision(10)<<Ncl<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep<<endl;
  
  phi[0]=prephi[0];
  phi[1]=prephi[1];
  dN1*=0.1;
  N=Ncl;

  while(ep<=1.0){
    prephi[0]=phi[0];
    prephi[1]=phi[1];
    ofs<< setprecision(10)<<N<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep<<endl;
    RK4(dphidN, N, phi, dN1);
    ep=pow(phi[1],2.0)/2.0/(pow((sqrt((1.0/6.0)*((pow(phi[1],2.0))+pow((mm*phi[0]),2.0)))),2.0));
  }
  Ncl=N-dN1;
  ep=pow(prephi[1],2.0)/2.0/(pow((sqrt((1.0/6.0)*((pow(prephi[1],2.0))+pow((mm*prephi[0]),2.0)))),2.0));

  cout<< setprecision(10)<<Ncl<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep<<endl;
  
  phi[0]=prephi[0];
  phi[1]=prephi[1];
  dN1*=0.1;
  N=Ncl;

  while(ep<=1.0){
    prephi[0]=phi[0];
    prephi[1]=phi[1];
    ofs<< setprecision(10)<<N<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep<<endl;
    RK4(dphidN, N, phi, dN1);
    ep=pow(phi[1],2.0)/2.0/(pow((sqrt((1.0/6.0)*((pow(phi[1],2.0))+pow((mm*phi[0]),2.0)))),2.0));
  }
  Ncl=N-dN1;
  ep=pow(prephi[1],2.0)/2.0/(pow((sqrt((1.0/6.0)*((pow(prephi[1],2.0))+pow((mm*prephi[0]),2.0)))),2.0));

  cout<< setprecision(10)<<Ncl<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep<<endl;
  
  phi[0]=prephi[0];
  phi[1]=prephi[1];
  dN1*=0.1;
  N=Ncl;

  while(ep<=1.0){
    prephi[0]=phi[0];
    prephi[1]=phi[1];
    ofs<< setprecision(10)<<N<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep<<endl;
    RK4(dphidN, N, phi, dN1);
    ep=pow(phi[1],2.0)/2.0/(pow((sqrt((1.0/6.0)*((pow(phi[1],2.0))+pow((mm*phi[0]),2.0)))),2.0));
  }
  Ncl=N-dN1;
  ep=pow(prephi[1],2.0)/2.0/(pow((sqrt((1.0/6.0)*((pow(prephi[1],2.0))+pow((mm*prephi[0]),2.0)))),2.0));

  cout<< setprecision(10)<<Ncl<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep<<endl;
  
  phi[0]=prephi[0];
  phi[1]=prephi[1];
  dN1*=0.1;
  N=Ncl;

  while(ep<=1.0){
    prephi[0]=phi[0];
    prephi[1]=phi[1];
    ofs<< setprecision(10)<<N<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep<<endl;
    RK4(dphidN, N, phi, dN1);
    ep=pow(phi[1],2.0)/2.0/(pow((sqrt((1.0/6.0)*((pow(phi[1],2.0))+pow((mm*phi[0]),2.0)))),2.0));
  }
  Ncl=N-dN1;
  ep=pow(prephi[1],2.0)/2.0/(pow((sqrt((1.0/6.0)*((pow(prephi[1],2.0))+pow((mm*prephi[0]),2.0)))),2.0));

  cout<< setprecision(10)<<Ncl<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep<<endl;
  
  phi[0]=prephi[0];
  phi[1]=prephi[1];
  dN1*=0.1;
  N=Ncl;

  while(ep<=1.0){
    prephi[0]=phi[0];
    prephi[1]=phi[1];
    ofs<< setprecision(10)<<N<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep<<endl;
    RK4(dphidN, N, phi, dN1);
    ep=pow(phi[1],2.0)/2.0/(pow((sqrt((1.0/6.0)*((pow(phi[1],2.0))+pow((mm*phi[0]),2.0)))),2.0));
  }
  Ncl=N-dN1;
  ep=pow(prephi[1],2.0)/2.0/(pow((sqrt((1.0/6.0)*((pow(prephi[1],2.0))+pow((mm*prephi[0]),2.0)))),2.0));

  cout<< setprecision(10)<<Ncl<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep<<endl;
  
  phi[0]=prephi[0];
  phi[1]=prephi[1];
  dN1*=0.1;
  N=Ncl;

  while(ep<=1.0){
    prephi[0]=phi[0];
    prephi[1]=phi[1];
    ofs<< setprecision(10)<<N<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep<<endl;
    RK4(dphidN, N, phi, dN1);
    ep=pow(phi[1],2.0)/2.0/(pow((sqrt((1.0/6.0)*((pow(phi[1],2.0))+pow((mm*phi[0]),2.0)))),2.0));
  }
  Ncl=N-dN1;
  ep=pow(prephi[1],2.0)/2.0/(pow((sqrt((1.0/6.0)*((pow(prephi[1],2.0))+pow((mm*prephi[0]),2.0)))),2.0));

  cout<< setprecision(10)<<Ncl<<"   "<<prephi[0]<<"   "<<prephi[1]<<"  "<<ep<<endl;
  
  phi[0]=prephi[0];
  phi[1]=prephi[1];
  dN1*=0.1;
  N=Ncl;



  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
  // -------------------------------------
}

vector<double> dphidN(double N, vector<double> phi) {
  vector<double> dphidN(2);

  dphidN[0] = (phi[1]/(sqrt((1.0/6.0)*((pow(phi[1],2.0))+pow((mm*phi[0]),2.0)))));//dphi/dN=(pi/H)　ノイズなし
  dphidN[1] = -3.0*phi[1]-((pow(mm,2.0))*phi[0])/(sqrt((1.0/6.0)*((pow(phi[1],2.0))+pow((mm*phi[0]),2.0))));//dpi/dN=-3pi-(v'/H)

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