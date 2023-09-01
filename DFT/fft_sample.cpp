// g++ -std=c++11 -O2 -o fft_sample fft_sample.cpp

#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>
#include <sys/time.h>

using namespace std;

vector<complex<double>> dft(vector<complex<double>> signal);
vector<complex<double>> fft(vector<complex<double>> signal);

const string filename = "fft_sin.dat";
ofstream ofs(filename);

int main() {
    int Num = pow(2, 3);
    double knum = 2;
    vector<complex<double>> signal;
    for (int i = 0; i < Num; i++) {
        complex<double> number(cos(2*M_PI*i/Num*knum), 0);
        signal.push_back(number);
    }
    
    // ---------- start timer ----------
    struct timeval Nv;
    struct timezone Nz;
    double before, after;
    
    gettimeofday(&Nv, &Nz);
    before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
    // --------------------------------------

    // FFT
    vector<complex<double>> spectrumf = fft(signal);
    
    // ---------- stop timer ----------
    gettimeofday(&Nv, &Nz);
    after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
    cout << "FFT " << after - before << " sec." << endl;
    // -------------------------------------
   

    
    gettimeofday(&Nv, &Nz);
    before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;

    // DFT
    vector<complex<double>> spectrumd = dft(signal);
    
    gettimeofday(&Nv, &Nz);
    after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
    cout << "DFT " << after - before << " sec." << endl;


    // 離散フーリエ変換結果を複素数形式で表示
    for (int k = 0; k < signal.size(); ++k) {
        // cout << k << " " << (spectrumf[k]-spectrumd[k]).real() << " + " << (spectrumf[k]-spectrumd[k]).imag() << "i" << endl;
        // cout << "Frequency bin " << k << ": " << spectrumd[k].real() << " + " << spectrumd[k].imag() << "i" << endl;
        ofs << spectrumf[k].real() << " " << spectrumf[k].imag() << endl;
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