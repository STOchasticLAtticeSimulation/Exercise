// g++ -std=c++11 -O2 -o fft_sample fft_sample.cpp

#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>
#include <sys/time.h>
#include <random>

using namespace std;

const complex<double> II(0,1);

vector<complex<double>> dft(vector<complex<double>> signal);
vector<vector<complex<double>>> dft(vector<vector<complex<double>>> signal);
vector<vector<vector<complex<double>>>> dft(vector<vector<vector<complex<double>>>> signal);
vector<complex<double>> fft(vector<complex<double>> signal);
vector<vector<complex<double>>> fft(vector<vector<complex<double>>> signal);
vector<vector<vector<complex<double>>>> fft(vector<vector<vector<complex<double>>>> signal);

bool realpoint(int nx, int ny, int nz, int Num); // judge real point
bool complexpoint(int nx, int ny, int nz, int Num); // judge independent complex point

// random distribution
std::random_device seed;
std::mt19937 engine(seed());
std::normal_distribution<> dist(0., 1.);

int main() {
    int Num = pow(2, 5);
    vector<vector<vector<complex<double>>>> signal(Num, vector<vector<complex<double>>>(Num, vector<complex<double>>(Num, 0)));
    int count = 0;
    
    for (int i = 0; i < Num; i++) {
      for (int j = 0; j < Num; j++) {
	for (int k = 0; k < Num; k++) {
	  if (realpoint(i,j,k,Num)) {
	    signal[i][j][k] = dist(engine);
	    count++;
	  } else if (complexpoint(i,j,k,Num)) {
	    signal[i][j][k] = (dist(engine) + II*dist(engine))/sqrt(2);
	    count += 2;
	  }
	}
      }
    }

    std::cout << Num*Num*Num << ' ' << count << std::endl;


    /*
      
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
    */
}


bool realpoint(int nx, int ny, int nz, int Num) {
  return (nx==0||nx==Num/2) && (ny==0||ny==Num/2) && (nz==0||nz==Num/2);
}

bool complexpoint(int nx, int ny, int nz, int Num) {
  int nxt, nyt, nzt; // shifted index
  if (nx<=Num/2) {
    nxt = nx;
  } else {
    nxt = nx-Num;
  }

  if (ny<=Num/2) {
    nyt = ny;
  } else {
    nyt = ny-Num;
  }

  if (nz<=Num/2) {
    nzt = nz;
  } else {
    nzt = nz-Num;
  }

  return (1<=nxt && nxt!=Num/2 && nyt!=Num/2 && nzt!=Num/2) ||
    (nxt==Num/2 && nyt!=Num/2 && nzt!=Num/2) || (nxt!=Num/2 && nyt==Num/2 && nzt!=Num/2) || (nxt!=Num/2 && nyt!=Num/2 && nzt==Num/2) ||
    (nxt==0 && nyt!=Num/2 && 1<=nzt && nzt!=Num/2) ||
    (nxt==Num/2 && nyt==Num/2 && 1<=nzt && nzt!=Num/2) || (nxt==Num/2 && 1<=nyt && nyt!=Num/2 && nzt==Num/2) || (1<=nxt && nxt!=Num/2 && nyt!=Num/2 && nzt!=Num/2) ||
    ((nxt==0 || nxt==Num/2) && 1<=nyt && nyt!=Num/2 && (nzt==0 || nzt==Num/2));
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



