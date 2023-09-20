// g++ -std=c++11 -O2 -o fft_sample fft_sample.cpp

#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>
#include <sys/time.h>

using namespace std;

const complex<double> II(0,1);

vector<complex<double>> dft(vector<complex<double>> signal);
vector<vector<complex<double>>> dft(vector<vector<complex<double>>> signal);
vector<vector<vector<complex<double>>>> dft(vector<vector<vector<complex<double>>>> signal);
vector<complex<double>> fft(vector<complex<double>> signal);
vector<vector<complex<double>>> fft(vector<vector<complex<double>>> signal);
vector<vector<vector<complex<double>>>> fft(vector<vector<vector<complex<double>>>> signal);

int main() {
  int Num = pow(2, 8);
  double kx1 = 10;
  double kx2 = 29;
  double ky1 = 21;
  double ky2 = 1;
  double kz1 = 3;
  double kz2 = 18;
  vector<vector<vector<complex<double>>>> signal(Num, vector<vector<complex<double>>>(Num, vector<complex<double>>(Num, 0)));
  for (int i = 0; i < Num; i++) {
    for (int j = 0; j < Num; j++) {
      for (int k = 0; k < Num; k++) {
	signal[i][j][k] = 1./Num/Num/Num*(exp(2*M_PI*(kx1*i+ky1*j+kz1*k)/Num*II) + exp(2*M_PI*(kx2*i+ky2*j+kz2*k)/Num*II));
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
  vector<vector<vector<complex<double>>>> spectrumf = fft(signal);
  
  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  cout << "FFT " << after - before << " sec." << endl;
  // -------------------------------------
  
  
  
  gettimeofday(&Nv, &Nz);
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  
  // DFT
  vector<vector<vector<complex<double>>>> spectrumd = dft(signal);
  
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  cout << "DFT " << after - before << " sec." << endl;
  
  
  for (int i = 0; i < signal.size(); i++) {
    for (int j = 0; j < signal[0].size(); j++) {
      for (int k = 0; k < signal[0][0].size(); k++) {
	if (abs(spectrumf[i][j][k]) > 1e-10) {
	  cout << i << ' ' << j << ' ' << k << ' ' << abs(spectrumf[i][j][k]) << endl;
	}
      }
    }
  }
  
  return 0;
}


// DFT
vector<complex<double>> dft(vector<complex<double>> signal) {
    int N = signal.size();
    vector<complex<double>> spectrum(N);

    for (int k = 0; k < N; k++) {
        complex<double> sum(0.0, 0.0);
        for (int n = 0; n < N; ++n) {
            double angle = -2*M_PI*k*n/N;
            sum += signal[n] * polar(1.0, angle);
        }
        spectrum[k] = sum;
    }

    return spectrum;
}

vector<vector<complex<double>>> dft(vector<vector<complex<double>>> signal) {
  // dft index y
  for (vector<complex<double>> &e : signal) {
    e = dft(e);
  }

  // transpose
  vector<vector<complex<double>>> tmp = signal;
  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
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



