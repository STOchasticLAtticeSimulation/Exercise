#include "STOLAS.hpp"
#include "RK4.hpp"

STOLAS::STOLAS(std::string Model, double NF, std::string noisedir, int noisefileNo, double Bias, double NBias, double DNbias) {
  model = Model;
  Nf = NF;
  bias = Bias;
  Nbias = NBias;
  dNbias = DNbias;

  std::cout << "model : " << model << std::endl;
  
  noisefile.open(noisedir+noisefilename);
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

  std::cout << "noise data imported" << std::endl;
}
