#include "../source/STOLAS.hpp"

const std::string model = "chaotic";
const double Nf = 5.;
const std::string noisedir = "../source/noisedata";
const int noisefileNo = 0.;
const double bias = 10.;
const double Nbias = 3.;
const double dNbias = 0.1;

const double mm = 0.01;

double STOLAS::VV(double phi) {
  return mm*mm*phi*phi/2.;
}

double STOLAS::Vp(double phi) {
  return mm*mm*phi;
}


int main()
{
  STOLAS stolas(model,Nf,noisedir,noisefileNo,bias,Nbias,dNbias);
}
