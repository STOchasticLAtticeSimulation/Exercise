#include "../source/STOLAS_Sta.hpp"
#include <sys/time.h>

const std::string model = "Starobinsky";
const std::string sourcedir = "../source";
const double dN = 0.01;
const double sigma = 0.1;

const double H0 = 1e-5;
const double calPRIR = 8.5e-10;
const double Lambda = 4574; //1e+4;
const double Ap = sqrt(9./4/M_PI/M_PI*H0*H0*H0*H0*H0*H0/calPRIR);
const double Am = Ap/Lambda;
const double V0 = 3*H0*H0;
const std::vector<double> phii{0.0302,-5.45e-7};
const double phif = -0.0182;
const double bias = 10.;
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

double STOLAS::calPphi(double &N, std::vector<double> &phi, double N0, bool broken) {
  if (!broken) {
    return pow(hubble(phi[0],phi[1])/2./M_PI,2);
  } else {
    double alpha = exp(N-N0);
    return pow(hubble(phi[0],phi[1])/2./M_PI,2) *
      ((1 + pow(sigma,2))*(9*pow(1 + pow(alpha,2)*pow(sigma,2),2) - 
				18*Lambda*pow(1 + pow(alpha,2)*pow(sigma,2),2) + 
				pow(Lambda,2)*(9 + 18*pow(alpha,2)*pow(sigma,2) + 9*pow(alpha,4)*pow(sigma,4) + 
					       2*pow(alpha,6)*pow(sigma,6))) + 
	    3*(-3*(1 + (-1 + 4*alpha)*pow(sigma,2) - (-4 + alpha)*pow(alpha,3)*pow(sigma,4) + 
		   pow(alpha,4)*pow(sigma,6)) + Lambda*
	       (6 + 6*(-1 + 4*alpha)*pow(sigma,2) + 2*(14 - 5*alpha)*pow(alpha,3)*pow(sigma,4) + 
		2*(5 - 2*alpha)*pow(alpha,4)*pow(sigma,6)) + 
	       pow(Lambda,2)*(-3 + (3 - 12*alpha)*pow(sigma,2) + pow(alpha,3)*(-16 + 7*alpha)*pow(sigma,4) + 
			      pow(alpha,4)*(-7 + 4*alpha)*pow(sigma,6)))*cos(2*(-1 + alpha)*sigma) + 
	    6*sigma*(-3*pow(-1 + Lambda,2) + pow(alpha,4)*(3 - 10*Lambda + 7*pow(Lambda,2))*pow(sigma,4) - 
		     3*alpha*pow(-1 + Lambda,2)*(-1 + pow(sigma,2)) - 
		     pow(alpha,3)*(3 - 7*Lambda + 4*pow(Lambda,2))*pow(sigma,2)*(-1 + pow(sigma,2)) + 
		     pow(alpha,5)*(-1 + Lambda)*Lambda*pow(sigma,4)*(-1 + pow(sigma,2)))*sin(2*sigma - 2*alpha*sigma))/
      (2.*pow(alpha,6)*pow(Lambda,2)*pow(sigma,6));
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
  stolas.compaction();
  //stolas.animation();
  //stolas.powerspec();

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  std::cout << std::endl;
  // -------------------------------------
}
