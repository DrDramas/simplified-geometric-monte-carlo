CXXFLAGS := -O2 -pthread $(shell root-config --cflags)
LDLIBS   := $(shell root-config --libs)

all: GetParticleDistributions SGMC

GetParticleDistributions: GetParticleDistributions.cpp
	$(CXX) $(CXXFLAGS) $< $(LDLIBS) -o $@

SGMC: SGMC.cpp
	$(CXX) $(CXXFLAGS) $< $(LDLIBS) -o $@

clean:
	rm -f GetParticleDistributions SGMC