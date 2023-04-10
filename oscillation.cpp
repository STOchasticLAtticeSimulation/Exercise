// 単振動を Runge--Kutta 4 で解いて t, x, xdot を "oscillation.dat" に出力する

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

void RK4forHO(double &t, double &x, double &xdot, double dt); // t, x, xdot を渡すと単振動 EoM に従い dt だけ t, x, xdot を更新する

// ------------ パラメータ ----------------- //
const string filename = "oscillation.dat"; // 出力ファイル名
const double tf = 10; // 終了時刻
const double dt = 0.01; // 時間刻み
// ----------------------------------------- //

// 変数の初期値
const double x0 = 1.0;
const double xdot0 = 0.0;

int main()
{
    ofstream ofs("oscillator.dat");
    double t = 0.0;
    double x = x0;
    double xdot = xdot0;
    while (t<tf) {
        cout << t << ' ' << x << ' ' << xdot << endl;
        ofs << t << ' ' << x << ' ' << xdot << endl;
        RK4forHO(t, x, xdot, dt);
    }
}

void RK4forHO(double &t, double &x, double &xdot, double dt){
    vector<double> kx(4);
    vector<double> kxdot(4);
    double a[4][4],b[4],c[4];

    for(int i=0;i<=3;i++){
        kx[i]=0.0;
        kxdot[i]=0.0;
        for(int j=0;j<=3;j++){
            a[i][j]=0;
        }
    }

    a[1][0]=1./2;  a[2][1]=1./2;  a[3][2]=1.;
    b[0]=1./6;     b[1]=1./3;     b[2]=1./3;    b[3]=1./6; 
    c[0]=0;        c[1]=1./2;     c[2]=1./2;    c[3]=1;

    for(int i=0;i<=3;i++){
        vector<double> X(4);
        vector<double> XDOT(4);
        
        for(int j=0;j<=3;j++){
            X[i]+=dt*a[i][j]*kx[j];
            XDOT[i]+=dt*a[i][j]*kxdot[j];

            kx[i]=xdot+XDOT[i];
            kxdot[i]=-(x+X[i]);
        }
    }
    x+=dt*(b[0]*kx[0]+b[1]*kx[1]+b[2]*kx[2]+b[3]*kx[3]);
    xdot+=dt*(b[0]*kxdot[0]+b[1]*kxdot[1]+b[2]*kxdot[2]+b[3]*kxdot[3]);
    
    t+=dt;
}
