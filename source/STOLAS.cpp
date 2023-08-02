#include "STOLAS.hpp"
#include "RK4.hpp"

STOLAS::STOLAS(std::string Model, double NF, std::string noisedir, int noisefileNo, double Bias, double NBias, double DNbias) {

#ifdef _OPENMP
  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;
#endif
  
  model = Model;
  Nf = NF;
  bias = Bias;
  Nbias = NBias;
  dNbias = DNbias;

  std::cout << "model : " << model << std::endl;
  
  noisefile.open(noisedir + std::string("/") + noisefilename);
  std::string str;
  std::stringstream ss;
  double dd;
  while (std::getline(noisefile, str)) {
    std::vector<double> vv;
    ss << str;
    while (!ss.eof()) {
      ss >> dd;
      vv.push_back(dd);
    }
    noisedata.push_back(vv);
  }

  NL = cbrt(noisedata.size());
  dN = Nf / (noisedata[0].size()-2);
  std::cout << "Noise data imported. Box size is " << NL << ". Time step is " << dN << " e-folds. Simulation ends at " << Nf << " e-folds." << std::endl;
}


/*
void STOLAS::dNmap() {
  
}
*/


double STOLAS::ep(double phi, double pi) {
  double HH = hubble(phi,pi);
  return pi*pi/2./HH/HH;
}

double STOLAS::hubble(double phi, double pi) {
  return sqrt((pi*pi/2. + VV(phi))/3.);
}

std::vector<double> STOLAS::dphidN(double N, std::vector<double> phi) {
  std::vector<double> dphidN(2);

  double xx = phi[0]; // phi
  double pp = phi[1]; // pi
  double HH = hubble(xx,pp);

  dphidN[0] = pp/HH;
  dphidN[1] = -3*pp - Vp(xx)/HH;
  
  return dphidN;
}

std::vector<double> STOLAS::dphidNbias(double N, std::vector<double> phi, std::vector<double> pos) {
  std::vector<double> dphidN(2);

  double xx = phi[0]; // phi
  double pp = phi[1]; // pi
  double HH = hubble(xx,pp);

  // bias
  double GaussianFactor = 1./dNbias /sqrt(2.*M_PI) * exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);
  double GaussianBias = bias * GaussianFactor;

  int ix = (int)(pos[0]+0.5) / (NL*NL);
  int iy = (int)(pos[0]+0.5) / NL;
  int iz = (int)(pos[0]+0.5) % NL;

  double r = sqrt(ix*ix + iy*iy + iz*iz);
  double ksigma = 2.*M_PI*sigma*exp(N)/NL;

  if (r!=0) {
    GaussianBias *= sin(ksigma*r)/ksigma/r;
  }

  dphidN[0] = pp/HH + HH/2./M_PI * GaussianBias;
  dphidN[1] = -3*pp - Vp(xx)/HH;
  
  return dphidN;
}

void STOLAS::RK4M(std::function<std::vector<double>(double, std::vector<double>)> dphidN, double dw, double &N, std::vector<double> &phi, double dN) {
  double HH = hubble(phi[0],phi[1]);
  
  RK4<std::vector<double>>(dphidN,N,phi,dN);
  phi[0] += HH/2./M_PI * dw;
}

void STOLAS::RK4Mbias(std::function<std::vector<double>(double, std::vector<double>, std::vector<double>)> dphidNbias, double dw, double &N, std::vector<double> &phi, double dN, std::vector<double> pos) {
  double HH = hubble(phi[0],phi[1]);

  RK4<std::vector<double>>(dphidNbias,N,phi,dN,pos);
  phi[0] += HH/2./M_PI * dw;
}
