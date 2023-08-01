#ifndef INCLUDED_STOLAS_

#define INCLUDED_STOLAS_

#include <cmath>
#include <functional>
#include <random>
#include "vec_op.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

// ------------ USER DEFINE ----------------
double VV(double phi); // potential V
double Vp(double phi); // \partial_\phi V
// -----------------------------------------

// ------------ User may change --------------------------
const double sigma = 0.1; // coarse-graining param. 
const double kdx = 0.1; // ksigma Deltax / 2pi.
// -------------------------------------------------------

void RK4(std::function<std::vector<double>(double, std::vector<double>)> dphidN, double &N, std::vector<double> &phi, double dN); // update e-folds &N and &phi = {phi, pi} by time step dN following EoM dphidN in Runge--Kutta 4
std::vector<double> dphidN(double N, std::vector<double> phi); // EoM
double ep(double phi, double pi); // epsilonH = -Hdot / H^2
double hubble(double phi, double pi); // Hubble param.

double Ncl(std::vector<double> phi, double N, double dN, double Nprec); // classical efolds Ncl from i.c. phi & N to the end of inflation epsilonH = 1 at the precision Nprec with the initial time step dN

template <class T>
void EulerM(std::function<T(double, const T&, double)> dNlist, std::function<T(double, const T&, double)> dwlist, double &N, T &x, double dN); // Euler--Maruyama
template <class T>
void RK4M(std::function<T(double, const T&, double)> dNlistRK4, std::function<T(double, const T&, double)> dwlist, double &N, T &x, double dN);

std::vector<std::vector<std::vector<std::vector<double>>>> dNlist(double N, std::vector<std::vector<std::vector<std::vector<double>>>> xif, double dN); // coeff. of dN
std::vector<std::vector<std::vector<std::vector<double>>>> dNlistRK4(double N, std::vector<std::vector<std::vector<std::vector<double>>>> xif, double dN);
std::vector<std::vector<std::vector<std::vector<double>>>> dwlist(double N, std::vector<std::vector<std::vector<std::vector<double>>>> xif, double dN); // coeff. of dW

// random distribution
std::random_device seed;
std::mt19937 engine(seed());
std::normal_distribution<> dist(0., 1.);


void RK4(std::function<std::vector<double>(double, std::vector<double>)> dphidN, double &N, std::vector<double> &phi, double dN) {
  std::vector<double> kx[4]; // 4-stage slopes kx
  double a[4][4],b[4],c[4]; // Butcher

  // -------------- initialise kx, a, b, c --------------- //
  for(int i=0;i<=3;i++){
    kx[i].assign(phi.size(), 0.);
    
    for(int j=0;j<=3;j++){
      a[i][j]=0.;
    }
  }
  
  a[1][0]=1./2;  a[2][1]=1./2;  a[3][2]=1.;
  b[0]=1./6;     b[1]=1./3;     b[2]=1./3;    b[3]=1./6; 
  c[0]=0;        c[1]=1./2;     c[2]=1./2;    c[3]=1;
  // ----------------------------------------------------- //
  

  std::vector<double> X = phi; // position at i-stage
    
  for(int i=0;i<=3;i++){
    X = phi; // initialise X
    
    for(int j=0;j<=3;j++){
      for (size_t xi = 0, size = phi.size(); xi < size; xi++) {
	X[xi] += dN * a[i][j] * kx[j][xi];
      }

      kx[i] = dphidN(N+c[i]*dN, X); // EoM
    }
  }

  for (size_t xi = 0, size = phi.size(); xi < size; xi++) {
    phi[xi] += dN*(b[0]*kx[0][xi] + b[1]*kx[1][xi] + b[2]*kx[2][xi] + b[3]*kx[3][xi]);
  }
  
  N+=dN;
}

std::vector<double> dphidN(double N, std::vector<double> phi) {
  std::vector<double> dphidN(2);

  double xx = phi[0]; // phi
  double pp = phi[1]; // pi
  double HH = hubble(xx,pp);

  dphidN[0] = pp/HH;
  dphidN[1] = -3*pp - Vp(xx)/HH;
  
  return dphidN;
}

double ep(double phi, double pi) {
  double HH = hubble(phi,pi);
  return pi*pi/2./HH/HH;
}

double hubble(double phi, double pi) {
  return sqrt((pi*pi/2. + VV(phi))/3.);
}

double Ncl(std::vector<double> phi, double N, double dN, double Nprec){
  double dN1 = dN;
  std::vector<double> prephi(2);

  while(dN1 >= Nprec) {
    while(ep(phi[0],phi[1])<=1.0){
      prephi[0]=phi[0];
      prephi[1]=phi[1];
      RK4(dphidN, N, phi, dN1);
    }
    N -= dN1;
    
    phi[0]=prephi[0];
    phi[1]=prephi[1];
    dN1*=0.1;
  }
  return N;
}


template <class T>
void EulerM(std::function<T(double, const T&, double)> dNlist, std::function<T(double, const T&, double)> dwlist, double &N, T &x, double dN)
{
  T xem = x;
  x += dNlist(N, xem, dN);
  x += dwlist(N, xem, dN);
  N += dN;
}

template <class T>
void RK4M(std::function<T(double, const T&, double)> dNlistRK4, std::function<T(double, const T&, double)> dwlist, double &N, T &x, double dN)
{
  T xem = x;
  x = dNlistRK4(N, xem, dN);
  x += dwlist(N, xem, dN);
  N += dN;
}


std::vector<std::vector<std::vector<std::vector<double>>>> dNlist(double N, std::vector<std::vector<std::vector<std::vector<double>>>> xif, double dN) {
  std::vector<std::vector<std::vector<std::vector<double>>>> dNlist = xif;
  
  for(size_t i = 0; i < xif.size(); i++) for(size_t j = 0; j < xif[i].size(); j++) for(size_t k = 0; k < xif[i][j].size(); k++) {
	std::vector<double> phi(2);
	phi[0] = xif[i][j][k][0];
	phi[1] = xif[i][j][k][1];
	dNlist[i][j][k][0] = dphidN(N, phi)[0] * dN;
	dNlist[i][j][k][1] = dphidN(N, phi)[1] * dN;
      }
  
  return dNlist;
}

std::vector<std::vector<std::vector<std::vector<double>>>> dNlistRK4(double N, std::vector<std::vector<std::vector<std::vector<double>>>> xif, double dN) {
  std::vector<std::vector<std::vector<std::vector<double>>>> dNlistRK4 = xif;
  
  for(size_t i = 0; i < xif.size(); i++) for(size_t j = 0; j < xif[i].size(); j++) for(size_t k = 0; k < xif[i][j].size(); k++) {
	std::vector<double> phi(2);
	phi[0] = xif[i][j][k][0];
	phi[1] = xif[i][j][k][1];
	RK4(dphidN, N, phi, dN);
	N -= dN;
	dNlistRK4[i][j][k][0] = phi[0];
	dNlistRK4[i][j][k][1] = phi[1];
      }

  return dNlistRK4;
}

std::vector<std::vector<std::vector<std::vector<double>>>> dwlist(double N, std::vector<std::vector<std::vector<std::vector<double>>>> xif, double dN) {
  std::vector<std::vector<std::vector<std::vector<double>>>> dwlist = xif; 

  // make correlation noise
  int NL = xif.size();
  double ksigma = 2. * M_PI * sigma * exp(N) / NL;
  double dtheta = 2. * M_PI * kdx / NL / ksigma;
  int divth = int(M_PI / dtheta);
  std::vector<double> divph(divth);
  std::vector<double> thetai(divth);
  std::vector<double> dphi(divth);
  std::vector<double> dOmegai(divth);

  // ------------ make noise ----------------
  std::vector<std::vector<double>> Omegalist;
  //double bias = 10;
  //double dNGaussianInv = 1./0.05;
  for (int n = 0; n < divth; n++) {
    thetai[n] = (n + 0.5) * dtheta;
    dphi[n] = 2. * M_PI * kdx / NL / ksigma / sin(thetai[n]);
    divph[n] = int(2. * M_PI / dphi[n]);
    for (int l = 0; l < divph[n]; l++){
      double dOmegai = sin(thetai[n]) * dtheta * dphi[n];
      //double GaussianFactor = dNGaussianInv / sqrt(2.*M_PI) * exp(-0.5*(N-Nbias)*(N-Nbias)*dNGaussianInv*dNGaussianInv) * sqrt(dN*dOmegai/M_PI) * 0.5;
      //double GaussianBais = bias * GaussianFactor;
      //normal_distribution<> dist1(GaussianBais, 1.);
      double phii = l * dphi[n];
      //Omegalist.push_back({thetai[n], dphi[n], phii, sqrt(dN) * dist1(engine)});
      Omegalist.push_back({thetai[n], dphi[n], phii, sqrt(dN) * dist(engine)});
    } 
  }
  // -----------------------------------------
  
  // -------------- add noise ----------------
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0; i < xif.size(); i++) for(size_t j = 0; j < xif[i].size(); j++) for(size_t k = 0; k < xif[i][j].size(); k++) {
	dwlist[i][j][k][0] = 0.; // initialize
	dwlist[i][j][k][1] = 0.;
	for (size_t n = 0; n < Omegalist.size(); n++){
	  double thet = Omegalist[n][0];
	  double dpht = Omegalist[n][1];
	  double phit = Omegalist[n][2];
	  double dwt = Omegalist[n][3];
	  double dOmegai = sin(thet) * dtheta * dpht;
	  //double ksx = (i+0.5) * sin(thet)*cos(phit) + (j+0.5) * sin(thet)*sin(phit) + (k+0.5) * cos(thet);
	  double ksx = i * sin(thet)*cos(phit) + j * sin(thet)*sin(phit) + k * cos(thet);
	  ksx *= ksigma;
	  dwlist[i][j][k][0] += 0.5 * sqrt(dOmegai / M_PI) * (cos(ksx) - sin(ksx)) * dwt;
	}
	dwlist[i][j][k][0] *= 0.5 * hubble(xif[i][j][k][0], xif[i][j][k][1]) / M_PI; 
      }
  // -----------------------------------------

  return dwlist;
}




#endif
