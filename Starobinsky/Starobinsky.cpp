#include "../source/STOLAS_Sta.hpp"
#include <sys/time.h>

const std::string model = "Starobinsky";
const std::string sourcedir = "../source";
const double dN = 0.01;

const double H0 = 1e-5;
const double calPRIR = 8.5e-10;
const double Lambda = 1e+4;
const double Ap = sqrt(9./4/M_PI/M_PI*H0*H0*H0*H0*H0*H0/calPRIR);
const double Am = Ap/Lambda;
const double V0 = 3*H0*H0;
const std::vector<double> phii{0.0302,-5.45e-7};
const double phif = -0.0182;
const double bias = 0.;
const double Nbias = 4.;
const double dNbias = 1.;


double STOLAS::VV(double phi) {
  if (phi > 0) {
    return V0 + Ap*phi;
  } else {
    return V0 + Am*phi;
  }
}

double STOLAS::Vp(double phi) {
  if (phi > 0) {
    return Ap;
  } else {
    return Am;
  }
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
  
  STOLAS stolas(model,dN,sourcedir,noisefileNo,phii,bias,Nbias,dNbias,phif);

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

  if (stolas.Nfilefail()//||stolas.cmpfilefail()
      ) {
    std::cout << "The export file couldn't be opened. 'mkdir data'" << std::endl;
    return -1;
  }

  
  stolas.dNmap();
  //stolas.compaction();
  //stolas.animation();
  stolas.powerspec();

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  // -------------------------------------
}
