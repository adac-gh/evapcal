TARGETS= evapcal

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

CXXFLAGS = -Wall -O2 $(ROOTFLAGS)
CXXLIBS = $(ROOTLIBS)

all:$(TARGETS)

evapcal:evapcal.o
	g++ -o $@ evapcal.o $(CXXLIBS)

.C.o:
	g++ -c $(CXXFLAGS) $<

clean:
	rm *.o $(TARGETS)
