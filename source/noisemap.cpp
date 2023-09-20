#include <sys/time.h>
#include <iostream>

#include "fft.hpp"

int main()
{
  int Num = pow(2, 7);
  double kx1 = 10;
  double kx2 = 29;
  double ky1 = 21;
  double ky2 = 1;
  double kz1 = 3;
  double kz2 = 18;
  std::vector<std::vector<std::vector<std::complex<double>>>> signal(Num, std::vector<std::vector<std::complex<double>>>(Num, std::vector<std::complex<double>>(Num, 0)));
  for (int i = 0; i < Num; i++) {
    for (int j = 0; j < Num; j++) {
      for (int k = 0; k < Num; k++) {
	signal[i][j][k] = 1./Num/Num/Num*(std::polar(1., 2*M_PI*(kx1*i+ky1*j+kz1*k)/Num) + std::polar(1., 2*M_PI*(kx2*i+ky2*j+kz2*k)/Num));
      }
    }
  }
  
  // ---------- start timer ----------
  struct timeval Nv;
  struct timezone Nz;
  double before, after;
  
  gettimeofday(&Nv, &Nz);
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  // --------------------------------------


  // FFT
  std::vector<std::vector<std::vector<std::complex<double>>>> spectrumf = fft(signal);
  
  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << "FFT " << after - before << " sec." << std::endl;
  // -------------------------------------
  
  gettimeofday(&Nv, &Nz);
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;

  
  // DFT
  std::vector<std::vector<std::vector<std::complex<double>>>> spectrumd = dft(signal);
  
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << "DFT " << after - before << " sec." << std::endl;
  
  
  for (int i = 0; i < signal.size(); i++) {
    for (int j = 0; j < signal[0].size(); j++) {
      for (int k = 0; k < signal[0][0].size(); k++) {
	if (abs(spectrumf[i][j][k]) > 1e-10) {
	  std::cout << i << ' ' << j << ' ' << k << ' ' << abs(spectrumf[i][j][k]) << std::endl;
	}
      }
    }
  }
  
  return 0;
}
