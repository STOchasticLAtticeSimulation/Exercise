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
const double dN = 0.01;
const std::string sourcedir = "../source";
const std::vector<double> phii{6.08726, -4.43529e-6};//phii{6.82869,-2.85674e-6};
const double bias = 0.0; //7.5;
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



int main(int argc, char* argv[])
{
  if (argc!=2) {
    std::cout << "Specify the noise file number correctly." << std::endl;
    return -1;
  }
  
  // ---------- start timer ----------
  struct timeval Nv;
  struct timezone Nz;
  double before, after;
  
  gettimeofday(&Nv, &Nz);
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  // --------------------------------------

  int noisefileNo = atoi(argv[1]);
  
  STOLAS stolas(model,dN,sourcedir,noisefileNo,phii,bias,Nbias,dNbias);

  if (!stolas.checknoisefile()) {
    std::cout << "The noise file couldn't be opened." << std::endl;
    return -1;
  }

  if (!stolas.checkbiasfile()) {
    std::cout << "The bias file couldn't be opened." << std::endl;
    return -1;
  }

  if (!stolas.noisebiassize()) {
    std::cout << "The box sizes of the noise and the bias are inconsistent." << std::endl;
    return -1;
  }

  if (stolas.Nfilefail()||stolas.wfilefail()) {
    std::cout << "The export file couldn't be opened. 'mkdir data'" << std::endl;
    return -1;
  }

  /*
  if (stolas.Hfilefail()||stolas.pifilefail()) {
    std::cout << "Caution: export files for animation couldn't be opened." << std::endl;
  }
  */

  stolas.dNmap();
  stolas.powerspec();

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  // -------------------------------------
}
