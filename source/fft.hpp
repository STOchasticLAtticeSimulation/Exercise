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


#endif
