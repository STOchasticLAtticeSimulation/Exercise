#include "STOLAS.hpp"
#include "vec_op.hpp"


STOLAS::STOLAS(std::string Model, double NF, std::string noisedir, int noisefileNo, std::vector<double> Phii, double Bias, double NBias, double DNbias) {

#ifdef _OPENMP
  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;
#endif
  
  model = Model;
  Nf = NF;
  phii = Phii;
  Nfilename = Nfileprefix + std::to_string((int)Nf) + std::string("_noise_") + std::to_string(noisefileNo) + std::string(".dat");
  bias = Bias;
  Nbias = NBias;
  dNbias = DNbias;

  std::cout << "model : " << model << std::endl;
  
  noisefile.open(noisedir + std::string("/") + noisefilename);
  std::string str;
  double dd;
  while (std::getline(noisefile, str)) {
    std::vector<double> vv;
    std::stringstream ss(str);
    while (!ss.eof()) {
      ss >> dd;
      vv.push_back(dd);
    }
    vv.pop_back();
    noisedata.push_back(vv);
  }

  NL = cbrt(noisedata.size());
  dN = Nf / (noisedata[0].size()-1);
  std::cout << "Noise data imported. Box size is " << NL << ". Time step is " << dN << " e-folds. Simulation ends at " << Nf << " e-folds." << std::endl;
}


void STOLAS::dNmap() {
  std::ofstream Nfile(Nfilename);
  int complete = 0;
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<NL*NL*NL; i++) {
    double N=0;
    std::vector<double> phi = phii;
    for (size_t n=0; n<noisedata[i].size(); n++) {
      RK4Mbias(N,phi,dN,noisedata[i][n],i);
    }

    double dN1 = dN;
    std::vector<double> prephi(2);
    while (dN1 >= Nprec) {
      while (ep(phi[0],phi[1])<=1.) {
	prephi[0] = phi[0];
	prephi[1] = phi[1];
	RK4(N,phi,dN);
      }
      N -= dN1;

      phi[0] = prephi[0];
      phi[1] = prephi[1];
      dN1 *= 0.1;
    }

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      Nfile << i << ' ' << N << std::endl;
      complete++;
      std::cout << "\r" << complete << "/" << NL*NL*NL << std::flush;
    }
  }
  std::cout << std::endl;
}


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

std::vector<double> STOLAS::dphidNbias(double N, std::vector<double> phi, int pos) {
  std::vector<double> dphidNbias(2);

  double xx = phi[0]; // phi
  double pp = phi[1]; // pi
  double HH = hubble(xx,pp);

  // bias
  double GaussianFactor = 1./dNbias /sqrt(2.*M_PI) * exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);
  double GaussianBias = bias * GaussianFactor;

  int ix = pos / (NL*NL);
  int iy = pos / NL;
  int iz = pos % NL;

  double r = sqrt(ix*ix + iy*iy + iz*iz);
  double ksigma = 2.*M_PI*sigma*exp(N)/NL;

  if (r!=0) {
    GaussianBias *= sin(ksigma*r)/ksigma/r;
  }

  dphidNbias[0] = pp/HH + HH/2./M_PI * GaussianBias;
  dphidNbias[1] = -3*pp - Vp(xx)/HH;
  
  return dphidNbias;
}

void STOLAS::RK4(double &t, std::vector<double> &x, double dt) {
  std::vector<double> kx[4]; // 4-stage slopes kx
  double a[4][4],b[4],c[4]; // Butcher

  // -------------- initialise kx, a, b, c --------------- //
  for(int i=0;i<=3;i++){
    kx[i] = x;
    vec_op::init(kx[i]);
    
    for(int j=0;j<=3;j++){
      a[i][j]=0.;
    }
  }
  
  a[1][0]=1./2;  a[2][1]=1./2;  a[3][2]=1.;
  b[0]=1./6;     b[1]=1./3;     b[2]=1./3;    b[3]=1./6; 
  c[0]=0;        c[1]=1./2;     c[2]=1./2;    c[3]=1;
  // ----------------------------------------------------- //
  

  std::vector<double> X = x; // position at i-stage
    
  for(int i=0;i<=3;i++){
    X = x; // initialise X
    
    for(int j=0;j<=3;j++){
      X += dt * a[i][j] * kx[j];
      kx[i] = dphidN(t,X);
    }
  }

  t += dt;
  x += dt*(b[0]*kx[0] + b[1]*kx[1] + b[2]*kx[2] + b[3]*kx[3]);
}

void STOLAS::RK4bias(double &t, std::vector<double> &x, double dt, int pos) {
  std::vector<double> kx[4]; // 4-stage slopes kx
  double a[4][4],b[4],c[4]; // Butcher

  // -------------- initialise kx, a, b, c --------------- //
  for(int i=0;i<=3;i++){
    kx[i] = x;
    vec_op::init(kx[i]);
    
    for(int j=0;j<=3;j++){
      a[i][j]=0.;
    }
  }
  
  a[1][0]=1./2;  a[2][1]=1./2;  a[3][2]=1.;
  b[0]=1./6;     b[1]=1./3;     b[2]=1./3;    b[3]=1./6; 
  c[0]=0;        c[1]=1./2;     c[2]=1./2;    c[3]=1;
  // ----------------------------------------------------- //
  

  std::vector<double> X = x; // position at i-stage
    
  for(int i=0;i<=3;i++){
    X = x; // initialise X
    
    for(int j=0;j<=3;j++){
      X += dt * a[i][j] * kx[j];
      kx[i] = dphidNbias(t,X,pos);
    }
  }

  t += dt;
  x += dt*(b[0]*kx[0] + b[1]*kx[1] + b[2]*kx[2] + b[3]*kx[3]);
}

void STOLAS::RK4M(double &N, std::vector<double> &phi, double dN, double dw) {
  double HH = hubble(phi[0],phi[1]);
  
  RK4(N,phi,dN);
  phi[0] += HH/2./M_PI * dw;
}

void STOLAS::RK4Mbias(double &N, std::vector<double> &phi, double dN, double dw, int pos) {
  double HH = hubble(phi[0],phi[1]);

  RK4bias(N,phi,dN,pos);
  phi[0] += HH/2./M_PI * dw;
}
