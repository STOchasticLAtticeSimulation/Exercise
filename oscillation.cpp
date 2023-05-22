// 単振動を Runge--Kutta 4 で解いて t, x, xdot を "oscillation.dat" に出力する

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sys/time.h>
#include <functional>

using namespace std;

void RK4(function<vector<double>(double, vector<double>)> dxdt,
	 double &t, vector<double> &x, double dt); // t と (x, xdot) を渡すと EoM dxdt に従い dt だけ t, (x, xdot) を更新する
vector<double> dxdt(double t, vector<double> x); // t と (x, xdot) の関数としての EoM dxdt. 今回は単振動

// ------------ パラメータ ----------------- //
const string filename = "oscillation.dat"; // 出力ファイル名
const double tf = 10; // 終了時刻
const double dt = 0.01; // 時間刻み
// ----------------------------------------- //

// 変数の初期値
const double x0 = 0.0;
const double xdot0 = 1.0;

int main()
{
  // ---------- start timer ----------
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  // --------------------------------------
  
  
  ofstream ofs(filename);
  double t = 0.0;
  vector<double> x{x0, xdot0};
  while (t<tf) {
    cout << t << ' ' << x[0] << ' ' << x[1] << endl;
    ofs << t << ' ' << x[0] << ' ' << x[1] << endl;
    RK4(dxdt, t, x, dt);
  }


  // ---------- stop timer ----------
  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
  // -------------------------------------
}

vector<double> dxdt(double t, vector<double> x) {
  vector<double> dxdt(2);

  dxdt[0] = x[1];
  dxdt[1] = -x[0];

  return dxdt;
}

void RK4(function<vector<double>(double, vector<double>)> dxdt, double &t, vector<double> &x, double dt) {
  vector<double> kx[4]; // 4段分の勾配ベクトル kx
  double a[4][4],b[4],c[4]; // ブッチャー係数

  // -------------- kx, a, b, c の初期化 --------------- //
  for(int i=0;i<=3;i++){
    kx[i].assign(x.size(), 0.);
    
    for(int j=0;j<=3;j++){
      a[i][j]=0.;
    }
  }
  
  a[1][0]=1./2;  a[2][1]=1./2;  a[3][2]=1.;
  b[0]=1./6;     b[1]=1./3;     b[2]=1./3;    b[3]=1./6; 
  c[0]=0;        c[1]=1./2;     c[2]=1./2;    c[3]=1;
  // ----------------------------------------------------- //
  

  vector<double> X = x; // i段目の位置用のベクトル
    
  for(int i=0;i<=3;i++){
    X = x; // X の初期化
    
    for(int j=0;j<=3;j++){
      for (size_t xi = 0, size = x.size(); xi < size; xi++) {
	X[xi] += dt * a[i][j] * kx[j][xi];
      }

      kx[i] = dxdt(t+c[i]*dt, X); // EoM
    }
  }

  for (size_t xi = 0, size = x.size(); xi < size; xi++) {
    x[xi] += dt*(b[0]*kx[0][xi] + b[1]*kx[1][xi] + b[2]*kx[2][xi] + b[3]*kx[3][xi]);
  }
  
  t+=dt;
}
