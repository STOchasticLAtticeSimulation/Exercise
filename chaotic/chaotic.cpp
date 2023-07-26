#include "../source/STOLAS.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/time.h>

const double mm = 1e-2;
const std::string filename = "Ncl_chaotic.dat";
const std::string filename_c = "Lattice_chaotic.dat";
const std::string filename_f = "field_chaotic.dat";
const double dN = 0.01;
const int EXPSTEP = 10;
const double Nf = 5;
const double phi0 = 15.;
const double pi0 = -0.1*mm*mm;
const int NL = 9;
const double NPREC = 1e-7;

// useful macro
#define LOOP for(int i = 0; i < NL; i++) for(int j = 0; j < NL; j++) for(int k = 0; k < NL; k++)


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

    EulerM<std::vector<std::vector<std::vector<std::vector<double>>>>>(dNlist, dwlist, N, x, dN); // Euler-Maruyama 1step <vector<vector<vector<vector<double>>>>>
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
