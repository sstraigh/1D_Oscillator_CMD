CXX = g++
CXXFLAGS = -lgsl -lgslcblas

all: 1D-oscillator Position-Probability-Distribution Create-Position-Autocorrelation

Position-Probability-Distribution: create_pos_histogram.o
	$(CXX) -o Position-Probability-Distribution create_pos_histogram.o

Create-Position-Autocorrelation: centroid_pos_auto-tcf.o
	$(CXX) -o Create-Position-Autocorrelation centroid_pos_auto-tcf.o

centroid_pos_auto-tcf.o: centroid_pos_auto-tcf.cpp
	$(CXX) -c centroid_pos_auto-tcf.cpp

create_pos_histogram.o: create_pos_histogram.cpp
	$(CXX) -c create_pos_histogram.cpp

1D-oscillator: thermostat.o 1D_oscillator.o particle.o necklace.o transition_matrix.o gaussian.o
	$(CXX) $(CXXFLAGS) -o 1D_oscillator 1D_oscillator.o thermostat.o particle.o necklace.o transition_matrix.o gaussian.o

thermostat.o: thermostat.cpp
	$(CXX) -c thermostat.cpp

particle.o: particle.cpp thermostat.h
	$(CXX) -c particle.cpp

necklace.o: necklace.cpp particle.h thermostat.h gaussian.h
	$(CXX) -c necklace.cpp

transition_matrix.o: transition_matrix.cpp
	$(CXX) -c transition_matrix.cpp

gaussian.o: gaussian.cpp
	$(CXX) -c gaussian.cpp

1D-oscillator.o: 1D_oscillator.cpp thermostat.h particle.h necklace.h transition_matrix.h gaussian.h
	$(CXX) $(CXXFLAGS) -c 1D_oscillator.cpp

clean:
	rm -f 1D_oscillator Position-Probability-Distribution *.o *.gch
