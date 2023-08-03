MODEL = chaotic
SOURCE = ../source

CXX := g++-13
#CXX := g++

CXXFLAGS := -std=c++11 -O2 -fopenmp
#CXXFLAGS := -Xpreprcessor -fopenmp -lomp -std=c++11 -O2

$(MODEL): $(MODEL).o $(SOURCE)/STOLAS.o
	$(CXX) $(CXXFLAGS) -o $(MODEL) $(MODEL).o $(SOURCE)/STOLAS.o

$(MODEL).o: $(SOURCE)/STOLAS.hpp
STOLAS.o: $(SOURCE)/STOLAS.hpp $(SOURCE)/RK4.hpp $(SOURCE)/vec_op.hpp

clean:
	$(RM) *.o
	$(RM) $(SOURCE)/*.o
	$(RM) $(MODEL)
