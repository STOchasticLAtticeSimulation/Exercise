#include "../source/STOLAS.hpp"

const double mm = 1e-2;

double VV(double phi) {
  return mm*mm*phi*phi/2.;
}

double Vp(double phi) {
  return mm*mm*phi;
}

int main()
{}
