#ifndef INCLUDED_fft_hpp_
#define INCLUDED_fft_hpp_

#define _USR_MATH_DEFINES
#include <cmath>
#include <complex>
#include <vector>

const std::complex<double> II(0,1);

std::vector<std::complex<double>> dft(std::vector<std::complex<double>> signal); // 1-dim discrete Fourier trs.
std::vector<std::vector<std::complex<double>>> dft(std::vector<std::vector<std::complex<double>>> signal); // 2-dim discrete Fourier trs.
std::vector<std::vector<std::vector<std::complex<double>>>> dft(std::vector<std::vector<std::vector<std::complex<double>>>> signal); // 3-dim discrete Fouriere trs.
std::vector<std::complex<double>> fft(std::vector<std::complex<double>> signal); // 1-dim fast Fourier trs.
std::vector<std::vector<std::complex<double>>> fft(std::vector<std::vector<std::complex<double>>> signal); // 2-dim fast Fourier trs.
std::vector<std::vector<std::vector<std::complex<double>>>> fft(std::vector<std::vector<std::vector<std::complex<double>>>> signal); // 3-dim fast Fourier trs.


// DFT
std::vector<std::complex<double>> dft(std::vector<std::complex<double>> signal) {
  int N = signal.size();
  std::vector<std::complex<double>> spectrum(N);
  
  for (int k = 0; k < N; k++) {
    std::complex<double> sum(0.0, 0.0);
    for (int n = 0; n < N; ++n) {
      double angle = -2*M_PI*k*n/N;
      sum += signal[n] * polar(1.0, angle);
    }
    spectrum[k] = sum;
  }
  
  return spectrum;
}

std::vector<std::vector<std::complex<double>>> dft(std::vector<std::vector<std::complex<double>>> signal) {
  // dft index y
  for (std::vector<std::complex<double>> &e : signal) {
    e = dft(e);
  }

  // transpose
  std::vector<std::vector<std::complex<double>>> tmp = signal;
  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[x].size(); y++) {
      tmp[x][y] = signal[y][x];
    }
  }
  signal = tmp;

  // dft index x
  for (vector<complex<double>> &e : signal) {
    e = dft(e);
  }

  // transpose
  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
      tmp[x][y] = signal[y][x];
    }
  }

  return tmp;
}

vector<vector<vector<complex<double>>>> dft(vector<vector<vector<complex<double>>>> signal) {
  // dft indices y and z
  for (vector<vector<complex<double>>> &e : signal) {
    e = dft(e);
  }

  // (x,y,z) to (y,z,x)
  vector<vector<vector<complex<double>>>> tmp = signal;
  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
      for (size_t z = 0; z < signal[0][0].size(); z++) {
	tmp[y][z][x] = signal[x][y][z];
      }
    }
  }
  signal = tmp;
  
  // dft index x
  for (vector<vector<complex<double>>> &e1 : signal) {
    for (vector<complex<double>> &e2 : e1) {
      e2 = dft(e2);
    }
  }

  // (y,z,x) to (x,y,z)
  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
      for (size_t z = 0; z < signal[0][0].size(); z++) {
	tmp[x][y][z] = signal[y][z][x];
      }
    }
  }

  return tmp;
}


// Cooley-Tukey FFT
vector<complex<double>> fft(vector<complex<double>> signal) {
  int N = signal.size();
  if (N<2) return signal; // サイズが N=2 になるまで分割
  
  vector<complex<double>> even(N/2);
  vector<complex<double>> odd(N/2);
  for (int i=0; i<N/2; i++) {
    even[i] = signal[2*i];
    odd[i] = signal[2*i+1];
  }
  
  even = fft(even);
  odd = fft(odd);
  
  vector<complex<double>> spectrum(N);
  for (int k=0; k<N/2; k++) {
    complex<double> oddW = polar(1.0, -2*M_PI*k/N) * odd[k];
    spectrum[k] = even[k] + oddW;
    spectrum[k+N/2] = even[k] - oddW;
  }
  
  return spectrum;
}

vector<vector<complex<double>>> fft(vector<vector<complex<double>>> signal) {
  // fft index y
  for (vector<complex<double>> &e : signal) {
    e = fft(e);
  }

  // transpose
  vector<vector<complex<double>>> tmp = signal;
  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
      tmp[x][y] = signal[y][x];
    }
  }
  signal = tmp;

  // fft index x
  for (vector<complex<double>> &e : signal) {
    e = fft(e);
  }

  // transpose
  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
      tmp[x][y] = signal[y][x];
    }
  }

  return tmp;
}

vector<vector<vector<complex<double>>>> fft(vector<vector<vector<complex<double>>>> signal) {
  // fft indices y and z
  for (vector<vector<complex<double>>> &e : signal) {
    e = fft(e);
  }

  // (x,y,z) to (y,z,x)
  vector<vector<vector<complex<double>>>> tmp = signal;
  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
      for (size_t z = 0; z < signal[0][0].size(); z++) {
	tmp[y][z][x] = signal[x][y][z];
      }
    }
  }
  signal = tmp;
  
  // fft index x
  for (vector<vector<complex<double>>> &e1 : signal) {
    for (vector<complex<double>> &e2 : e1) {
      e2 = fft(e2);
    }
  }

  // (y,z,x) to (x,y,z)
  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
      for (size_t z = 0; z < signal[0][0].size(); z++) {
	tmp[x][y][z] = signal[y][z][x];
      }
    }
  }

  return tmp;
}

#endif
