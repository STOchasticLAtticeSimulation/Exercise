#include "../source/STOLAS.hpp"
#include <sys/time.h>

const double mm = 0.01;

const std::string model = "chaotic";
const double Nf = 5.5;
const std::string noisedir = "../source/noisedata";
const int noisefileNo = 0.;
const std::vector<double> phii{15.,-0.1*mm*mm};
const double bias = 10.;
const double Nbias = 2;
const double dNbias = 0.01;


double STOLAS::VV(double phi) {
  return mm*mm*phi*phi/2.;
}

double STOLAS::Vp(double phi) {
  return mm*mm*phi;
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
