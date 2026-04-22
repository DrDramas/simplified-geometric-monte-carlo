# Shared by every target
CXXFLAGS_COMMON := -O2 $(shell root-config --cflags)
LDLIBS          := $(shell root-config --libs)

# Parallel stages need pthreads
CXXFLAGS_PAR    := $(CXXFLAGS_COMMON) -pthread

all: RecoverBeamSpot GetParticleDistributions SGMC

RecoverBeamSpot: RecoverBeamSpot.cpp
	$(CXX) $(CXXFLAGS_COMMON) $< $(LDLIBS) -o $@

GetParticleDistributions: GetParticleDistributions.cpp
	$(CXX) $(CXXFLAGS_PAR) $< $(LDLIBS) -o $@

SGMC: SGMC.cpp
	$(CXX) $(CXXFLAGS_PAR) $< $(LDLIBS) -o $@

clean:
	rm -f RecoverBeamSpot GetParticleDistributions SGMC