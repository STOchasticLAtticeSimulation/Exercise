#include "../source/STOLAS.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/time.h>

const double AW = 0.02;
const double BW = 1;
const double CW = 0.04;
const double DW = 0;
const double GW = 3.076278e-2;
const double RW = 0.7071067;
const double CALV = 1000;
const double W0 = 12.35;
const double CUP = 0.0382;

const std::string filename = "Ncl_inflection.dat"; // 出力ファイル名(Ncl)
const std::string filename_c = "Lattice_inflection.dat"; // 出力ファイル名(曲率ゆらぎ)
const std::string filename_f = "field_inflection.dat"; // 出力ファイル名(phi, pi)
const double dN = 0.01;
const int EXPSTEP = 10;
const double Nf = 5;
const double phi0 = 3.60547;
const double pi0 = -2.37409e-7;
const int NL = 9;
const double NPREC = 1e-7;

// useful macro
#define LOOP for(int i = 0; i < NL; i++) for(int j = 0; j < NL; j++) for(int k = 0; k < NL; k++)

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
  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;
#endif

  // output
  int numsteps = 0;
  
  std::ofstream ofs(filename);
  std::ofstream ofs_c(filename_c);
  std::ofstream ofs_f(filename_f);
  ofs << std::setprecision(10);
  ofs_c << std::setprecision(10);
  ofs_f << std::setprecision(10);

  std::vector<double> xi{phi0,pi0};
  std::vector<std::vector<double>> x1(NL, xi); // dim 1
  std::vector<std::vector<std::vector<double>>> x2(NL, x1); // dim 2
  std::vector<double> v1(NL,0);
  std::vector<std::vector<double>> v2(NL,v1);
  std::vector<std::vector<std::vector<std::vector<double>>>> x(NL, x2); // dim 3

  // calculate average energy
  std::vector<std::vector<std::vector<double>>> de(NL, v2); // dim 3
  double average_e;

  double N = 0.;// e-foldings

  while (N < Nf) {
    // save the data
    if(numsteps % EXPSTEP == 0 && N < Nf){
      average_e = 0.; // initialize
      LOOP average_e += 3. * hubble(x[i][j][k][0], x[i][j][k][1]) * hubble(x[i][j][k][0], x[i][j][k][1]);
      average_e /= NL*NL*NL; // 平均エネルギー
      std::cout << N << ' ' << x[0][0][0][0] << ' ' << x[0][0][0][1] << ' ';
      std::vector<double> phi(2);
      LOOP{
        phi[0] = x[i][j][k][0];
        phi[1] = x[i][j][k][1];
        de[i][j][k] = 3. * hubble(phi[0], phi[1]) * hubble(phi[0], phi[1]) - average_e;
        ofs_c << de[i][j][k] / phi[1] / phi[1] / 3. << ' ';
      }
      std::cout << de[0][0][0]/phi[1]/phi[1]/3. << std::endl;
      ofs_c << std::endl;
    }
    numsteps++;

    RK4M<std::vector<std::vector<std::vector<std::vector<double>>>>>(dNlistRK4, dwlist, N, x, dN); // Euler-Maruyama 1step <vector<vector<vector<vector<double>>>>>
  }

  LOOP ofs_f << x[i][j][k][0] << ' ' << x[i][j][k][1] << std::endl; // 最終的な場の値を出力
  

  //------------- Ncl -------------
#ifdef _OPENMP
#pragma omp parallel for
#endif
  LOOP{
    std::vector<double> phi{x[i][j][k][0], x[i][j][k][1]};
    
    double NCL = Ncl(phi,N,dN,NPREC);
#ifdef _OPENMP
#pragma omp critical
#endif
    ofs << i+NL*j+NL*NL*k << ' ' << NCL << std::endl;
  }

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  // -------------------------------------
}
