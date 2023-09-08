#include "../source/STOLAS.hpp"
#include <sys/time.h>

const double AW = 0.02;
const double BW = 1;
const double CW = 0.04;
const double DW = 0;
const double GW = 3.076278e-2;
const double RW = 0.7071067;
const double CALV = 1000;
const double W0 = 12.35;
const double CUP = 0.0382;

const std::string model = "inflection";
const double Nf = 5.51;
const std::string noisedir = "../source/noisedata";
const int noisefileNo = 0.;
const std::vector<double> phii{3.60547,-2.37409e-7};
const double bias = 0.; //10.;
const double Nbias = 2;
const double dNbias = 0.01;


double STOLAS::VV(double phi) {
  return W0*W0/CALV/CALV/CALV * ( CUP/pow(CALV,1./3) + AW/(exp(phi/sqrt(3))-BW) - CW/exp(phi/sqrt(3))
				  + exp(2*phi/sqrt(3))/CALV*(DW-GW/(RW*exp(sqrt(3)*phi)/CALV+1)) );
}

double STOLAS::Vp(double phi) {
  return exp(-phi/sqrt(3)) * ( CW - AW*exp(2*phi/sqrt(3))/(BW-exp(phi/sqrt(3)))/(BW-exp(phi/sqrt(3)))
			       + 3*exp(2*sqrt(3)*phi)*GW*RW/(CALV+exp(sqrt(3)*phi)*RW)/(CALV+exp(sqrt(3)*phi)*RW)
			       + 2*exp(sqrt(3)*phi)*(DW-CALV*GW/(CALV+exp(sqrt(3)*phi)*RW))/CALV )
    * W0*W0/sqrt(3)/CALV/CALV/CALV;
}



int main()
{
  // ---------- start timer ----------
  struct timeval Nv;
  struct timezone Nz;
  double before, after;
  
  gettimeofday(&Nv, &Nz);
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  // --------------------------------------
  
  STOLAS stolas(model,Nf,noisedir,noisefileNo,phii,bias,Nbias,dNbias);

  stolas.dNmap();

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  // -------------------------------------
}
