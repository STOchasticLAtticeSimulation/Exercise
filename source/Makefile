CXX := g++

CXXFLAGS := -I/usr/local/opt/libomp/include -Xpreprocessor -fopenmp -std=c++11 -O2
LDFLAGS:= -L/usr/local/opt/libomp/lib

noisemap: fft.hpp
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o noisemap noisemap.cpp

clean:
	$(RM) noisemap
