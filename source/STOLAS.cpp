#include "STOLAS.hpp"
#include "vec_op.hpp"
#include "fft.hpp"

// useful macro
#define LOOP for(int i = 0; i < NL; i++) for(int j = 0; j < NL; j++) for(int k = 0; k < NL; k++)


STOLAS::STOLAS(std::string Model, double DN, std::string sourcedir, int noisefileNo, std::vector<double> Phii, double Bias, double NBias, double DNbias) {

#ifdef _OPENMP
  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;
#endif
  
  model = Model;
  dN = DN;
  phii = Phii;
  bias = Bias;
  Nbias = NBias;
  dNbias = DNbias;

  noisefile.open(sourcedir + std::string("/") + noisefilename + std::to_string(noisefileNo) + std::string(".dat"));
  noisefilefail = noisefile.fail();
  biasfile.open(sourcedir + std::string("/") + biasfilename);
  biasfilefail = biasfile.fail();
  
  if (!noisefile.fail() && !biasfile.fail() && !Nfile.fail()) {
    std::cout << "model : " << model << std::endl;
    
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

    while (std::getline(biasfile, str)) {
      std::vector<double> vv;
      std::stringstream ss(str);
      while (!ss.eof()) {
	ss >> dd;
	vv.push_back(dd);
      }
      vv.pop_back();
      biasdata.push_back(vv);
    }

    NL = cbrt(noisedata.size());
    std::cout << "Noise/Bias data imported. Box size is " << NL << "." << std::endl;
    Nfile.open(Nfileprefix + std::to_string(NL) + std::string("_") + std::to_string(noisefileNo) + std::string(".dat"));
    Hfile.open(Hfileprefix + std::to_string(NL) + std::string("_") + std::to_string(noisefileNo) + std::string(".dat"));
    pifile.open(pifileprefix + std::to_string(NL) + std::string("_") + std::to_string(noisefileNo) + std::string(".dat"));
    wfile.open(wfileprefix + std::to_string(NL) + std::string("_") + std::to_string(noisefileNo) + std::string(".dat"));
    powfile.open(powfileprefix + std::to_string(NL) + std::string("_") + std::to_string(noisefileNo) + std::string(".dat"));
    Nmap3D = std::vector<std::vector<std::vector<std::complex<double>>>>(NL, std::vector<std::vector<std::complex<double>>>(NL, std::vector<std::complex<double>>(NL, 0)));
  }
}


bool STOLAS::checknoisefile() {
  return !noisefilefail;
}

bool STOLAS::checkbiasfile() {
  return !biasfilefail;
}

bool STOLAS::Nfilefail() {
  return Nfile.fail();
}

bool STOLAS::Hfilefail() {
  return Hfile.fail();
}

bool STOLAS::pifilefail() {
  return pifile.fail();
}

bool STOLAS::wfilefail() {
  return wfile.fail();
}

bool STOLAS::powfilefail() {
  return powfile.fail();
}

void STOLAS::dNmap() {
  Nfile << std::setprecision(10);
  Hfile << std::setprecision(14);
  pifile << std::setprecision(14);
  wfile << std::setprecision(10);
  int complete = 0;

  std::vector<std::vector<double>> Hdata(noisedata[0].size(), std::vector<double>(NL*NL*NL,0));
  std::vector<std::vector<double>> pidata(noisedata[0].size(), std::vector<double>(NL*NL*NL,0));
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<NL*NL*NL; i++) {
    double N=0;
    std::vector<double> phi = phii;
    for (size_t n=0; n<noisedata[i].size(); n++) {
      Hdata[n][i] = pow(hubble(phi[0],phi[1]),2);
      pidata[n][i] = phi[1]*phi[1];
      //RK4Mbias(N,phi,dN,noisedata[i][n],i);
      RK4Mbias(N,phi,dN,noisedata[i][n],biasdata[i][n]);
    }

    double dN1 = dN;
    std::vector<double> prephi(2);
    while (dN1 >= Nprec) {
      while (ep(phi[0],phi[1])<=1.) {
	prephi[0] = phi[0];
	prephi[1] = phi[1];
	RK4(N,phi,dN1);
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
      std::cout << "\rLatticeSimulation : " << std::setw(3) << 100*complete/NL/NL/NL << "%" << std::flush;
    }

    //power spectrum
    int x=i/NL/NL ,y=(i%(NL*NL))/NL, z=i%NL;
    const std::complex<double> II(0,1);

    Nmap3D[x][y][z]=N;

  }
  std::cout << std::endl;

  for (size_t n=0; n<Hdata.size(); n++) {
    for (size_t i=0; i<Hdata[n].size(); i++) {
      Hfile << Hdata[n][i] << ' ';
      pifile << pidata[n][i] << ' ';
    }
    Hfile << std::endl;
    pifile << std::endl;
    std::cout << "\rAnimeDataExporting : " << std::setw(3) << 100*n/Hdata.size() << "%" << std::flush;
  }
  std::cout << "\rAnimeDataExporting : 100%" << std::endl;

  //calculation of weight
  double logw=0.;
  for(size_t n=0; n<noisedata[0].size(); n++){
    double N=n*dN;
    double Bias=bias *1./dNbias/sqrt(2*M_PI) * exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);
    logw+=-Bias*noisedata[0][n]*sqrt(dN)-(Bias*Bias*dN)/2;
  }
  wfile << logw << std::endl;

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

void STOLAS::powerspec(){
  powfile << std::setprecision(10);
  std::vector<std::vector<std::vector<std::complex<double>>>> Nk=fft(Nmap3D);
  //Nk = std::vector<std::vector<std::vector<double>>>(NL,std::vector<std::vector<double>>(NL,std::vector<double>(NL,0)));

  //std::vector<std::vector<std::vector<std::complex<double>>>> Nmap3Dfft = fft(Nmap3D);
  LOOP{
      int nxt, nyt, nzt; // shifted index
  if (i<=NL/2) {
    nxt = i;
  } else {
    nxt = i-NL;
  }

  if (j<=NL/2) {
    nyt = j;
  } else {
    nyt = j-NL;
  }

  if (k<=NL/2) {
    nzt = k;
  } else {
    nzt = k-NL;
  }
    
    double rk=nxt*nxt+nyt+nyt+nzt+nzt;
    powfile<< sqrt(rk) <<"     "<< norm(Nk[i][j][k]) << std::endl;
  }

}

/*
std::vector<double> STOLAS::dphidNbias(double N, std::vector<double> phi, int pos) {
  std::vector<double> dphidNbias(2);

  double xx = phi[0]; // phi
  double pp = phi[1]; // pi
  double HH = hubble(xx,pp);

  // bias
  double GaussianFactor = 1./dNbias /sqrt(2.*M_PI) * exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);
  double GaussianBias = bias * GaussianFactor;

  int ix = pos / (NL*NL);
  int iy = (pos%(NL*NL)) / NL;
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
*/

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

/*
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
  phi[0] += HH/2./M_PI * dw * sqrt(dN);
}
*/

void STOLAS::RK4Mbias(double &N, std::vector<double> &phi, double dN, double dw, double Bias //int pos
		      ) {
  double HH = hubble(phi[0],phi[1]);

  //RK4bias(N,phi,dN,pos);
  RK4(N,phi,dN);
  phi[0] += HH/2./M_PI * dw * sqrt(dN);

  double GaussianFactor = 1./dNbias/sqrt(2*M_PI) * exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);
  phi[0] += HH/2./M_PI * bias * Bias * GaussianFactor * dN;
}
