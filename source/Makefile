#MODEL = noisemap
MODEL = biasmap

#CXX := g++-13
CXX := g++

#CXXFLAGS := -std=c++11 -O2 -fopenmp
CXXFLAGS := -std=c++11 -O2 -Xpreprocessor -fopenmp -lomp

LDFLAGS := -L/opt/homebrew/opt/libomp/lib
CPPFLAGS := -I/opt/homebrew/opt/libomp/include

$(MODEL): $(MODEL).o
	$(CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $(MODEL) $(MODEL).o 

$(MODEL).o: vec_op.hpp fft.hpp

clean:
	$(RM) *.o
	$(RM) $(MODEL)
