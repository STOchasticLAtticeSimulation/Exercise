CXX := g++

CXXFLAGS := -std=c++11 -O2 #-I/usr/local/opt/libomp/include -Xpreprocessor -fopenmp -lomp
LDFLAGS:= -L/usr/local/opt/libomp/lib

noisemap: fft.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o noisemap noisemap.cpp

clean:
	$(RM) noisemap
