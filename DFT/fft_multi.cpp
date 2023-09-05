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
vector<complex<double>> fft(vector<complex<double>> signal);
vector<vector<complex<double>>> fft(vector<vector<complex<double>>> signal);

const string filename = "fft_exp.dat";
ofstream ofs(filename);

int main() {
    int Num = pow(2, 8);
    double kx1 = 10;
    double kx2 = 52;
    double ky1 = 21;
    double ky2 = 1;
    vector<vector<complex<double>>> signal(Num, vector<complex<double>>(Num, 0));
    for (int i = 0; i < Num; i++) {
      for (int j = 0; j < Num; j++) {
	signal[i][j] = 1./Num/Num*(exp(2*M_PI*(kx1*i+ky1*j)/Num*II) + exp(2*M_PI*(kx2*i+ky2*j)/Num*II));
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
    vector<vector<complex<double>>> spectrumf = fft(signal);
    
    // ---------- stop timer ----------
    gettimeofday(&Nv, &Nz);
    after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
    cout << "FFT " << after - before << " sec." << endl;
    // -------------------------------------
   

    
    gettimeofday(&Nv, &Nz);
    before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;

    // DFT
    vector<vector<complex<double>>> spectrumd = dft(signal);
    
    gettimeofday(&Nv, &Nz);
    after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
    cout << "DFT " << after - before << " sec." << endl;


    for (int i = 0; i < signal.size(); i++) {
      for (int j = 0; j < signal[0].size(); j++) {
	ofs << abs(spectrumf[i][j]) << ' ';
	if (abs(spectrumf[i][j]) > 1e-10) {
	  cout << i << ' ' << j << ' ' << abs(spectrumf[i][j]) << endl;
	}
      }
      ofs << endl;
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
  for (vector<complex<double>> &e : signal) {
    e = dft(e);
  }

  vector<vector<complex<double>>> tmp = signal;
  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
      tmp[x][y] = signal[y][x];
    }
  }
  signal = tmp;

  for (vector<complex<double>> &e : signal) {
    e = dft(e);
  }

  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
      tmp[x][y] = signal[y][x];
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
  for (vector<complex<double>> &e : signal) {
    e = fft(e);
  }

  vector<vector<complex<double>>> tmp = signal;
  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
      tmp[x][y] = signal[y][x];
    }
  }
  signal = tmp;

  for (vector<complex<double>> &e : signal) {
    e = fft(e);
  }

  for (size_t x = 0; x < signal.size(); x++) {
    for (size_t y = 0; y < signal[0].size(); y++) {
      tmp[x][y] = signal[y][x];
    }
  }

  return tmp;
}


